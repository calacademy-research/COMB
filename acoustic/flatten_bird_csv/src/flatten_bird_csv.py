"""Beam pipeline that converts logits CSV to Avro for BigQuery import.

Directory structure convention:

base_dir/
  input/
    birds.txt
    [**]/*.csv
  output/
    avro-*

Usage:

python3 flatten_bird_csv.py --base_dir=

Logits CSV has columns:

filename_stem,start_seconds,end_seconds,logit0,logit1,...

where start_second and end_seconds are the endpoints of a context window of
audio relative to the start of the file and where the ordering of logits
corresponds to the order of lines in birds.txt.

birds.txt is a tab-separated values file with no header line and the species
code in the second column.

The output consists of multiple Avro files (shards), where each record is a
prediction for a single context window. The fields of each record are: key,
bird, begin_seconds, logit.

To run on Google Cloud Dataflow, base_dir should be a gs:// path, --gcp_project=
and --gcp_region= should have values, and the GOOGLE_APPLICATION_CREDENTIALS
environment variable should be set to the path of service account JSON with
permissions to the storage bucket and the compute engine service account.

Dataflow does not support all Python versions, and dependencies can be tricky.
Best practice is:

python3.8 -m venv dataflow_env
. dataflow_env/bin/activate
pip install apache-beam[gcp]
"""

import logging
import os
from typing import Any, Dict, Iterable, Sequence

from absl import app
from absl import flags
from apache_beam.options import pipeline_options
import apache_beam as beam

FLAGS = flags.FLAGS

flags.DEFINE_string('base_dir', None, 'Parent directory of input/ and output/.')
flags.DEFINE_string('gcp_project', None,
                    'Cloud project name; enables DataflowRunner.')
flags.DEFINE_string('gcp_region', None,
                    'Cloud region. Required if gcp_project is set.')

AVRO_SCHEMA = {
    'namespace':
        'comb',
    'type':
        'record',
    'name':
        'Prediction',
    'fields': [
        {
            'name': 'key',
            'type': 'string'
        },
        {
            'name': 'bird',
            'type': 'string'
        },
        {
            'name': 'begin_seconds',
            'type': 'float'
        },
        {
            'name': 'logit',
            'type': 'float'
        },
    ]
}


def birds_txt_species_code(line: str) -> str:
  """Selects the lowercase species code from a line of birds.txt."""
  parts = line.strip().split('\t')
  code_parts = parts[1].split('"')
  assert len(code_parts[0]) == 0
  assert len(code_parts[2]) == 0
  return code_parts[1]


def flatten(line: str, bird_list: Iterable[str]) -> Iterable[Dict[str, Any]]:
  """Converts a predictions CSV line to per-species records.

  Args:
    line: A line from predictions CSV, which has the format (filename, begin,
      end, logit0, logit1, ...)

  Yields:
    A dict with the {key, bird, begin_seconds, logit} for each logit on this
    line.
  """
  parts = line.split(',')
  key = parts[0]
  begin_seconds = float(parts[1])
  logits = parts[3:]
  for code, logit in zip(bird_list, logits):
    logit = float(logit)
    yield {
        'key': key,
        'bird': code,
        'begin_seconds': begin_seconds,
        'logit': logit,
    }


def run(base_dir: str, options: pipeline_options.PipelineOptions) -> None:
  """Runs the CSV-to-Avro pipeline.

  Args:
    base_dir: Path to a directory containing input/ and output/ subdirectories.
      input/birds.txt should have tab-separated values with species codes in the
      second field. input/ should have any nested directory structure with logit
      CSV files ending in .csv.
    options: Runtime pipeline configuration.

  Returns:
    None
  """
  birds_txt = os.path.join(base_dir, 'input', 'birds.txt')
  csv_pattern = os.path.join(base_dir, 'input', '**', '*.csv')
  output = os.path.join(base_dir, 'output', 'avro')

  with beam.Pipeline(options=options) as pipeline:
    bird_list = beam.pvalue.AsSingleton(
        pipeline | ('ReadBirds' >> beam.io.textio.ReadFromText(birds_txt))
        | beam.Map(birds_txt_species_code) | beam.combiners.ToList())

    csv_lines = (
        pipeline | ('ReadCsv' >> beam.io.textio.ReadFromText(csv_pattern)))
    outputs = csv_lines | beam.ParDo(flatten, bird_list)
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
      options=options,
  )


if __name__ == '__main__':
  logging.getLogger().setLevel(logging.INFO)
  flags.mark_flag_as_required('base_dir')
  app.run(main)
