--
-- Template for database cleanup script.
-- Removes database tables corresponding to a particular visit.
--

-- association pipeline inputs
DROP TABLE IF EXISTS DiaSources_visit%(visitId)d;
DROP TABLE IF EXISTS MopsPreds_visit%(visitId)d;

-- association pipeline outputs
DROP TABLE IF EXISTS DiaSourceToObjectMatches_visit%(visitId)d;
DROP TABLE IF EXISTS MopsPredToDiaSourceMatches_visit%(visitId)d;
DROP TABLE IF EXISTS NewObjects_visit%(visitId)d;

-- temp tables created by the association pipeline store stage
DROP TABLE IF EXISTS BestMatch_visit%(visitId)d;
DROP TABLE IF EXISTS MatchedObjects_visit%(visitId)d;
DROP TABLE IF EXISTS NewObjects_visit%(visitId)d;

