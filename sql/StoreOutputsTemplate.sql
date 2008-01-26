--
-- Template for script which updates the DIASource and Object tables using 
-- the association pipeline outputs for a visit.
--
-- There are obvious opportunities for parallelization here. Exploiting them
-- is deferred to DC3, unless performance proves to be catastrophic.
--

-- Set objectId of each difference source to the id of the closest matching object
CREATE TABLE BestMatch_visit%(visitId)d LIKE InMemoryMatchPairTemplate;
INSERT INTO BestMatch_visit%(visitId)d
    SELECT first, second, MIN(distance)
    FROM DiaSourceToObjectMatches_visit%(visitId)d
    GROUP BY first;
UPDATE DiaSources_visit%(visitId)d AS s, BestMatch_visit%(visitId)d AS m
    SET s.objectId = m.second
    WHERE s.diaSourceId = m.first;
DROP TABLE BestMatch_visit%(visitId)d;

-- Set objectId of each difference source used to create an object
UPDATE DiaSources_visit%(visitId)d AS s, NewObjectIdPairs_visit%(visitId)d AS n
    SET s.objectId = n.second
    WHERE s.diaSourceId = n.first;

-- Finally, append difference sources to the historical DiaSource table
INSERT INTO %(diaSourceTable)s SELECT * FROM DiaSources_visit%(visitId)d;

-- Create list of matched object ids without duplicates
CREATE TABLE MatchedObjects_visit%(visitId)d LIKE InMemoryIdTemplate;
INSERT INTO MatchedObjects_visit%(visitId)d
    SELECT DISTINCT(second) FROM DiaSourceToObjectMatches_visit%(visitId)d;

-- Update latest observation time and observation count for objects with matches
UPDATE %(varObjectTable)s AS o, MatchedObjects_visit%(visitId)d AS m
    SET   o.latestObsTime = '%(visitTime)s',
          o.%(filterName)cNumObs = o.%(filterName)cNumObs + 1
    WHERE o.objectId = m.id;
UPDATE %(nonVarObjectTable)s AS o, MatchedObjects_visit%(visitId)d AS m
    SET   o.latestObsTime = '%(visitTime)s',
          o.%(filterName)cNumObs = o.%(filterName)cNumObs + 1
    WHERE o.objectId = m.id;
DROP TABLE MatchedObjects_visit%(visitId)d;

-- Create new objects from difference sources in an in-memory table. This provides both
-- a debugging aid and allows different filters to be handled with a single script.
CREATE TABLE NewObjects_visit%(visitId)d LIKE InMemoryObjectTemplate;

-- Insert records (all filter-specific fields set to NULL)
INSERT INTO NewObjects_visit%(visitId)d
    SELECT n.second, # objectId
           NULL,
           s.ra, # ra
           s.decl, # decl
           NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           '%(visitTime)s', # earliestObsTime
           '%(visitTime)s', # latestObsTime
           NULL, NULL, NULL, NULL, NULL,
           s.cx, # cx
           NULL,  
           s.cy, # cy
           NULL,  
           s.cz, # cz
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
           NULL, NULL, NULL, NULL
    FROM NewObjectIdPairs_visit%(visitId)d AS n
         JOIN DiaSources_visit%(visitId)d AS s ON (s.diaSourceId = n.first);

-- Set filter specific fields in the in-memory temp table
UPDATE NewObjects_visit%(visitId)d       AS o,
       NewObjectIdPairs_visit%(visitId)d AS n,
       DiaSources_visit%(visitId)d       AS s
    SET -- o.%(filterName)cMag          = NULL,
        -- o.%(filterName)cMagErr       = NULL,
        o.%(filterName)cErrA         = IF(s.raErr4detection > s.decErr4detection, s.raErr4detection, s.decErr4detection)/1000.0,
        o.%(filterName)cErrB         = IF(s.raErr4detection > s.decErr4detection, s.decErr4detection, s.raErr4detection)/1000.0,
        o.%(filterName)cErrTheta     = IF(s.raErr4detection > s.decErr4detection, 90, 0),
        o.%(filterName)cNumObs       = 1,
        o.%(filterName)cVarProb      = 100,
        -- o.%(filterName)cAmplitude    = NULL,
        -- o.%(filterName)cPeriod       = NULL,
        -- o.%(filterName)cPetroMag     = NULL,
        -- o.%(filterName)cPetroMagErr  = NULL,
        -- o.%(filterName)cIxx          = NULL,
        -- o.%(filterName)cIxxErr       = NULL,
        -- o.%(filterName)cIyy          = NULL,
        -- o.%(filterName)cIyyErr       = NULL,
        -- o.%(filterName)cIxy          = NULL,
        -- o.%(filterName)cIxyErr       = NULL,
        -- o.%(filterName)cTimescale    = NULL,
        o.%(filterName)cApMag        = s.apMag,
        o.%(filterName)cApMagErr     = s.apMagErr
        -- o.%(filterName)cIsoAreaImage = NULL,
        -- o.%(filterName)cMuMax        = NULL,
        -- o.%(filterName)cFluxRadius   = NULL,
        -- o.%(filterName)cFlags        = NULL
    WHERE o.objectId = n.second AND n.first = s.diaSourceId;

-- Finally, append the in-memory temp table to the variable object partition
INSERT INTO %(varObjectTable)s SELECT * FROM NewObjects_visit%(visitId)d;

DROP TABLE NewObjects_visit%(visitId)d;

