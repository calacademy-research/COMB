-- yARU
--
-- This query produces detection counts for each (year, species, point, visit).
--
-- An acoustic "visit" is defined to be the recording of one file
-- (GROUP BY predictions.key).
--
-- The visits are identified by their position in the sorted list of files for
-- the same year and point.
--
-- The result columns include detection counts at each of an arithmetic
-- sequence of logit thresholds.

SELECT
  year,
  bird_code,
  point_number,
  RANK() OVER (PARTITION BY year, bird_code, point_number ORDER BY file_start) AS visit_id,
  detection_count_m1_0,
  detection_count_m0_5,
  detection_count_0_0,
  detection_count_0_5,
  detection_count_1_0,
  detection_count_1_5,
  detection_count_2_0,
  detection_count_2_5,
  detection_count_3_0
FROM (
  SELECT

  -- id.vars
    bird_codes.code_4 AS bird_code,
    point.point_number AS point_number,
    EXTRACT(YEAR FROM begin_timestamp_millis) AS year,
    MIN(begin_timestamp_millis) AS file_start,

  -- measure.vars
  SUM(CAST(l.value > -1.0 AS INT64)) AS detection_count_m1_0,
  SUM(CAST(l.value > -0.5 AS INT64)) AS detection_count_m0_5,
  SUM(CAST(l.value > 0.0 AS INT64)) AS detection_count_0_0,
  SUM(CAST(l.value > 0.5 AS INT64)) AS detection_count_0_5,
  SUM(CAST(l.value > 1.0 AS INT64)) AS detection_count_1_0,
  SUM(CAST(l.value > 1.5 AS INT64)) AS detection_count_1_5,
  SUM(CAST(l.value > 2.0 AS INT64)) AS detection_count_2_0,
  SUM(CAST(l.value > 2.5 AS INT64)) AS detection_count_2_5,
  SUM(CAST(l.value > 3.0 AS INT64)) AS detection_count_3_0

  FROM comb.predictions,
  UNNEST(logits) AS l
  JOIN comb.bird_codes ON bird_codes.code_4 = l.key
  WHERE point.point_number > 0
  AND l.key NOT IN ('unknown', 'human', 'nonbird')

  GROUP BY bird_code, point_number, year, predictions.key
);
