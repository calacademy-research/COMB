"""Beam pipeline that converts logits CSV to Avro for BigQuery import.

Directory structure convention:

base_dir/  (passed in via --base_dir)
  input/
    birds.txt
    aru2point.csv
    latlong.csv
    [**]/*.csv
  output/
    avro-*

The files [**]/*.csv are logits CSV, with columns:

filename_stem,start_seconds,end_seconds,logit0,logit1,...

where start_second and end_seconds are the endpoints of a context window of
audio relative to the start of the file and where the ordering of logits
corresponds to the order of lines in birds.txt.

birds.txt is a tab-separated values file with no header line and the species
code in the second column.

aru2point.csv is a CSV file with no header and columns filename,point_number.
The point numbers must be integers. In cases where they are unknown, the line
should be omitted.

latlong.csv is a CSV file with no header and columns
point_number,latitude,longitude. The coordinates are in decimal degrees.

The output consists of multiple Avro files (shards), where each record is a
prediction for a single context window. See AVRO_SCHEMA for the schema.

To run on Google Cloud Dataflow, base_dir should be a gs:// path, --gcp_project=
and --gcp_region= should have values, and the GOOGLE_APPLICATION_CREDENTIALS
environment variable should be set to the path of service account JSON with
permissions to the Storage bucket and the Compute Engine service account.

Dataflow does not support all Python versions, and dependencies can be tricky.
Best practice is:

python3.8 -m venv dataflow_env
. dataflow_env/bin/activate
pip install apache-beam[gcp]
"""

import csv
import datetime
import io
import logging
import os
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

from absl import app
from absl import flags
from apache_beam.io import fileio
from apache_beam.options import pipeline_options
import apache_beam as beam

FLAGS = flags.FLAGS

flags.DEFINE_string('base_dir', None, 'Parent directory of input/ and output/.')
flags.DEFINE_string('gcp_project', None,
                    'Cloud project name; enables DataflowRunner.')
flags.DEFINE_string('gcp_region', None,
                    'Cloud region. Required if gcp_project is set.')
flags.DEFINE_float(
    'min_logit', None,
    'If set, logits less than this value will be dropped to sparsify output.')
flags.DEFINE_list(
    'codes_to_exclude', [],
    'All logits for species codes in this comma-separated list will be dropped.'
)

AVRO_SCHEMA = {
    'namespace':
        'comb',
    'type':
        'record',
    'name':
        'Predictions',
    'fields': [
        {
            'name': 'key',
            'type': 'string'
        },
        {
            'name': 'aru',
            'type': 'string'
        },
        {
            'name':
                'point',
            'type': [{
                'name':
                    'Point',
                'type':
                    'record',
                'fields': [
                    {
                        'name': 'point_number',
                        'type': 'int'
                    },
                    {
                        'name': 'latitude',
                        'type': 'float'
                    },
                    {
                        'name': 'longitude',
                        'type': 'float'
                    },
                ],
            }, 'null'],
        },
        {
            'name': 'begin_rel_file',
            'type': 'float'
        },
        {
            'name': 'end_rel_file',
            'type': 'float'
        },
        {
            'name': 'begin_timestamp_millis',
            'type': {
                'type': 'int',
                'logicalType': 'timestamp-millis',
            },
        },
        {
            'name': 'logits',
            'type': {
                'type': 'map',
                'values': 'float',
            },
        },
    ]
}

LatLong = Tuple[float, float]


def normalize_key(k):
  """Makes logits CSV and aru2point.csv keys correspond."""
  return os.path.splitext(k)[0]


def read_species_codes(tsv_file: fileio.ReadableFile) -> Sequence[str]:
  """Reads the mapping from class indices to species code.

  Args:
    tsv_file: Object with a .read_utf8() method that returns the entire contents
      of the species codes file in the input/ directory.

  Returns:
    Species codes ordered as they were in the file.
  """
  with io.StringIO(tsv_file.read_utf8()) as csvfile:
    reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
    return [row[1] for row in reader]


def parse_aru_line(aru_line: str) -> Tuple[str, int]:
  """Returns (key, point_number) from a line of aru2point.csv."""
  parts = aru_line.split(',')
  return normalize_key(parts[0]), int(parts[1])


def parse_latlong_line(latlong_line: str) -> Tuple[int, LatLong]:
  """Returns (point_number, (lat, long)) from a line of latlong.csv."""
  parts = latlong_line.split(',')
  point = int(parts[0])
  latitude = float(parts[1])
  longitude = float(parts[2])
  return point, (latitude, longitude)


class LogitCsvLineToAvroRecordFn(beam.DoFn):
  """Beam DoFn for joining logit CSV to point metadata."""

  def __init__(self, min_logit: float, codes_to_exclude: Iterable[str] = None):
    super(LogitCsvLineToAvroRecordFn, self).__init__()
    if not codes_to_exclude:
      codes_to_exclude = set()

    self._min_logit = min_logit
    self._codes_to_exclude = codes_to_exclude

  def process(self, logit_csv_line: str, species_codes: Sequence[str],
              aru_to_point: Dict[str, int],
              point_to_latlong: Dict[int, LatLong]) -> Iterable[Dict[str, Any]]:
    """Converts a logits CSV line to an Avro record.

    This does lookup joins from aru_to_point and point_to_latlong and adds this
    point metadata to the structured output.

    The logits are converted form a dense CSV row to a (species => logit) map,
    optinally sparsifying as specified in the constructor arguments.

    Args:
      logit_csv_line: A line from predictions CSV, which has the format
        (filename, begin, end, logit0, logit1, ...). The order of the logits
        should correspond to the order of species_codes.
      species_codes: List of strings whose order matches the ordering of class
        indices in the model output. These become keys in a (species => logit)
        map in the structured output.
      aru_to_point: Side input with the mapping from filenames to point number.
      point_to_latlong: Side input with the mapping of point number to
        coordinates.

    Yields:
      A dict serializable according to AVRO_SCHEMA.
    """
    parts = logit_csv_line.split(',')
    key = normalize_key(parts[0])
    begin_rel_file = float(parts[1])
    end_rel_file = float(parts[1])
    logits = parts[3:]

    key_parts = key.split('_')
    aru = key_parts[0]

    file_start_utc = datetime.datetime.strptime('_'.join(key_parts[1:]),
                                                '%Y%m%d_%H%M%S')
    begin_timestamp_millis = int(
        1000 * (file_start_utc +
                datetime.timedelta(seconds=begin_rel_file)).timestamp())

    point_number = aru_to_point.get(key, None)
    if point_number:
      latitude, longitude = point_to_latlong[point_number]
      point = {
          'point_number': point_number,
          'latitude': latitude,
          'longitude': longitude,
      }
    else:
      point = None

    logits_map: Dict[str, float] = {}
    for code, logit in zip(species_codes, logits):
      if code in self._codes_to_exclude:
        continue
      logit = float(logit)
      if (not self._min_logit) or logit >= self._min_logit:
        logits_map[code] = logit

    if logits_map:
      yield {
          'key': key,
          'aru': aru,
          'begin_timestamp_millis': begin_timestamp_millis,
          'begin_rel_file': begin_rel_file,
          'end_rel_file': end_rel_file,
          'point': point,
          'logits': logits_map,
      }


def run(base_dir: str, min_logit: Optional[float],
        codes_to_exclude: Iterable[str],
        options: pipeline_options.PipelineOptions) -> None:
  """Runs the CSV-to-Avro pipeline.

  Args:
    base_dir: Path to a directory containing input/ fitting the description in
      the module docstring.
    min_logit: If set, logits lower than this will be omitted from the output.
    codes_to_exclude: Species codes whose logits should always be omitted from
      the output.
    options: Runtime pipeline configuration.

  Returns:
    None
  """
  birds_txt_path = os.path.join(base_dir, 'input', 'birds.txt')
  arus_csv_path = os.path.join(base_dir, 'input', 'aru2point.csv')
  latlong_csv_path = os.path.join(base_dir, 'input', 'latlong.csv')
  csv_pattern = os.path.join(base_dir, 'input', '**', '*.csv')
  output = os.path.join(base_dir, 'output', 'avro')

  with beam.Pipeline(options=options) as pipeline:
    # Since to order of the lines in birds.txt determines the class indices, it
    # must be read all at once rather than with ReadFromText, which produces a
    # PCollection with no guaranteed order.
    species_codes = beam.pvalue.AsSingleton(
        pipeline | ('MatchBirds' >> fileio.MatchFiles(birds_txt_path))
        | fileio.ReadMatches() | beam.Map(read_species_codes))

    aru_to_point = beam.pvalue.AsDict(pipeline | (
        'ReadArus' >> beam.io.textio.ReadFromText(arus_csv_path))
                                      | beam.Map(parse_aru_line))

    point_to_latlong = beam.pvalue.AsDict(pipeline | (
        'ReadLatLong' >> beam.io.textio.ReadFromText(latlong_csv_path))
                                          | beam.Map(parse_latlong_line))

    csv_lines = (
        pipeline | ('ReadCsv' >> beam.io.textio.ReadFromText(csv_pattern)))
    outputs = csv_lines | beam.ParDo(
        LogitCsvLineToAvroRecordFn(min_logit, codes_to_exclude),
        species_codes,
        aru_to_point,
        point_to_latlong,
    )
    _ = outputs | beam.io.avroio.WriteToAvro(output, AVRO_SCHEMA)


def main(argv: Sequence[str]) -> None:
  del argv

  base_dir = FLAGS.base_dir

  if FLAGS.gcp_project:
    options = beam.options.pipeline_options.PipelineOptions(
        runner='DataflowRunner',
        project=FLAGS.gcp_project,
        region=FLAGS.gcp_region,
        job_name='flatten-bird-csv',
        temp_location=os.path.join(base_dir, 'tmp'),
    )
  else:
    options = beam.options.pipeline_options.PipelineOptions(
        runner='DirectRunner',
        direct_num_workers=0,
    )

  run(
      base_dir=base_dir,
      min_logit=FLAGS.min_logit,
      codes_to_exclude=set(FLAGS.codes_to_exclude),
      options=options,
  )


if __name__ == '__main__':
  logging.getLogger().setLevel(logging.INFO)
  flags.mark_flag_as_required('base_dir')
  app.run(main)
