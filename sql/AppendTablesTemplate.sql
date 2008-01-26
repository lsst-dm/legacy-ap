--
-- Template for accumulating per-visit table data.
-- Appends per-visit tables to global accumulator tables
--

INSERT INTO MopsPreds
    SELECT *, %(visitId)d FROM MopsPreds_visit%(visitId)d;
INSERT INTO DiaSourceToObjectMatches
    SELECT *, %(visitId)d FROM DiaSourceToObjectMatches_visit%(visitId)d;
INSERT INTO MopsPredToDiaSourceMatches
    SELECT *, %(visitId)d FROM MopsPredToDiaSourceMatches_visit%(visitId)d;
INSERT INTO NewObjectIdPairs
    SELECT *, %(visitId)d FROM NewObjectIdPairs_visit%(visitId)d;

