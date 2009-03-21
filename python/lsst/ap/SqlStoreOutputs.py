"""
Parameterized SQL statement templates which update the DIASource
and Object tables using association pipeline outputs for a visit
"""

mysqlStatements = [
    # Set objectId of each difference source to the id of the closest matching object
    """CREATE TEMPORARY TABLE _tmp_v%(visitId)d_BestMatch LIKE _tmpl_InMemoryMatchPair""",
    """ALTER TABLE _tmp_v%(visitId)d_DIASourceToObjectMatches ADD INDEX (first, distance)""",
    """INSERT INTO _tmp_v%(visitId)d_BestMatch
       SELECT a.first, a.second, a.distance
       FROM _tmp_v%(visitId)d_DIASourceToObjectMatches AS a LEFT OUTER JOIN
            _tmp_v%(visitId)d_DIASourceToObjectMatches AS b ON
            a.first = b.first AND
            (b.distance < a.distance OR (b.distance = a.distance AND b.second < a.second))
       WHERE b.first IS NULL""",
    """UPDATE _tmp_v%(visitId)d_DIASource AS s, _tmp_v%(visitId)d_BestMatch AS m
       SET s.objectId = m.second
       WHERE s.diaSourceId = m.first""",
    # Set objectId of each difference source used to create an object
    """UPDATE _tmp_v%(visitId)d_DIASource AS s, _tmp_v%(visitId)d_DIASourceToNewObject AS n
       SET s.objectId = n.second
       WHERE s.diaSourceId = n.first""",
    # Set id of sibling difference source (measured on other exposure in visit)
    """UPDATE _tmp_v%(visitId)d_DIASource
       SET diaSourceToId = diaSourceId ^ (1 << 14)""",
    # Append difference sources to the historical DiaSource table
    """INSERT INTO %(diaSourceTable)s SELECT * FROM _tmp_v%(visitId)d_DIASource""",
    # Update latest observation time and observation count for objects with matches
    """UPDATE %(varObjectTable)s AS o, _tmp_v%(visitId)d_BestMatch AS m
       SET   o.latestObsT = '%(dateObs)s',
             o.%(filter)cNumObs = o.%(filter)cNumObs + 1
       WHERE o.objectId = m.second""",
    """UPDATE %(nonVarObjectTable)s AS o, _tmp_v%(visitId)d_BestMatch AS m
       SET   o.latestObsT = '%(dateObs)s',
             o.%(filter)cNumObs = o.%(filter)cNumObs + 1
       WHERE o.objectId = m.second""",
    # Create new objects from difference sources in an in-memory table.
    """CREATE TEMPORARY TABLE _tmp_v%(visitId)d_NewObject LIKE _tmpl_InMemoryObject""",
    # Insert records (all filter-specific fields set to NULL)
    """INSERT INTO _tmp_v%(visitId)d_NewObject
       SELECT n.second, # objectId
           s.ra, # ra
           s.decl, # decl
           '%(dateObs)s', # earliestObsT
           '%(dateObs)s', # latestObsT
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL,
           NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL
       FROM _tmp_v%(visitId)d_DIASourceToNewObject AS n
       JOIN _tmp_v%(visitId)d_DIASource AS s ON (s.diaSourceId = n.first)""",
    # Set filter specific fields in the in-memory temp table
    # FIXME: copying flux values into magnitudes is clearly bogus 
    """UPDATE _tmp_v%(visitId)d_NewObject            AS o,
              _tmp_v%(visitId)d_DIASourceToNewObject AS n,
              _tmp_v%(visitId)d_DIASource            AS s
       SET o.%(filter)cMag         = s.psfFlux,
           o.%(filter)cMagErr      = s.psfFluxErr,
           o.%(filter)cIxx         = s.Ixx,
           o.%(filter)cIyy         = s.Iyy,
           o.%(filter)cIxy         = s.Ixy,
           o.%(filter)cNumObs      = 1
        WHERE o.objectId = n.second AND n.first = s.diaSourceId""",
    # Append the in-memory temp table to the variable object catalog
    """INSERT INTO %(varObjectTable)s SELECT * FROM _tmp_v%(visitId)d_NewObject"""
]

sqlStatements = { 'mysql': mysqlStatements }

