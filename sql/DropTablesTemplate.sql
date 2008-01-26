--
-- Template for database cleanup script.
-- Removes database tables corresponding to a particular visit.
--

-- drop association pipeline inputs
DROP TABLE IF EXISTS DiaSources_visit%(visitId)d;
DROP TABLE IF EXISTS MopsPreds_visit%(visitId)d;

-- drop association pipeline outputs
DROP TABLE IF EXISTS DiaSourceToObjectMatches_visit%(visitId)d;
DROP TABLE IF EXISTS MopsPredToDiaSourceMatches_visit%(visitId)d;
DROP TABLE IF EXISTS NewObjectIdPairs_visit%(visitId)d;

