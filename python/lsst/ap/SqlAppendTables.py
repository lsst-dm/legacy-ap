"""
Parameterized SQL statement templates for appending the contents
of per-visit tables to global "accumulator" tables.
"""

mysqlStatements = [
    """INSERT INTO _mops_Prediction
       SELECT *, %(visitId)d FROM _tmp_v%(visitId)d_Preds""",
    """INSERT INTO _ap_DIASourceToObjectMatches
       SELECT *, %(visitId)d FROM _tmp_v%(visitId)d_DIASourceToObjectMatches""",
    """INSERT INTO _ap_PredToDIASourceMatches
       SELECT *, %(visitId)d FROM _tmp_v%(visitId)d_PredToDIASourceMatches""",
    """INSERT INTO _ap_DIASourceToNewObject
       SELECT *, %(visitId)d FROM _tmp_v%(visitId)d_DIASourceToNewObject""",
]

sqlStatements = { 'mysql': mysqlStatements }
