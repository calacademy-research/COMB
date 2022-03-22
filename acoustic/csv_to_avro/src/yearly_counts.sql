-- Query used to materialize yearly_counts
--
-- yearly_counts is the data source for the Data Studio map visualization.
--
-- The `predictions` table from which this query selects is an import of the
-- output of csv_to_avro.py, without --min_logit or --codes_to_exclude flags
-- set.

SELECT
    l.key AS bird,
    point.point_number AS point,
    ANY_VALUE(point.latitude) AS latitude,
    ANY_VALUE(point.longitude) AS longitude,
    EXTRACT(YEAR FROM begin_timestamp_millis) AS year,
    SUM(CAST(l.value > 0 AS INT64)) AS detection_count,
    MAX(l.value) AS max_logit
FROM `comb.predictions_ge_-1`,
UNNEST(logits) AS l
GROUP BY bird, point, year;
