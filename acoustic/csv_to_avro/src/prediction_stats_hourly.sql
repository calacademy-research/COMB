-- Hourly rollup of dense logits
--
-- This query computes statistics of the logits grouped by species code and hour
-- for subsequent queryies where a flat table with (species_code, score) columns
-- may be convenient.
--
-- The `predictions` table from which this query selects is an import of the
-- output of csv_to_avro.py, without --min_logit or --codes_to_exclude flags
-- set.

SELECT
    l.key AS bird,
    TIMESTAMP_TRUNC(begin_timestamp_millis, HOUR) AS hour_utc,
    SUM(CAST(l.value > 0 AS INT64)) AS detection_count,
    MAX(l.value) AS max_logit,
    AVG(l.value) AS mean_logit,
    AVG(1 / (1 + EXP(-l.value))) AS max_probability,
    AVG(1 / (1 + EXP(-l.value))) AS mean_probability
FROM comb.predictions,
UNNEST(logits) AS l
GROUP BY bird, hour_utc;
