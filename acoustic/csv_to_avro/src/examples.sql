-- Gallery of exploratory queries

-- Per-species quantiles of hourly detection counts
SELECT
    bird,
    APPROX_QUANTILES(detection_count, 10) hourly_count_quantiles
FROM comb.prediction_stats_hourly
GROUP BY 1 ORDER BY 1;

-- Per-species hourly detections counts sorted by count
SELECT
  bird,
  detection_count,
  hour_utc
FROM comb.prediction_stats_hourly
WHERE bird != 'unknown'
ORDER BY detection_count DESC LIMIT 1000;
