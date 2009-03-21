"""
Parameterized SQL statement templates for dropping
temporary tables used during visit processing
"""

mysqlStatements = [
    # Association pipeline inputs
    """DROP TABLE IF EXISTS _tmp_v%(visitId)d_Preds""",
    """DROP TABLE IF EXISTS _tmp_v%(visitId)d_DIASource""",
    # Association pipeline outputs
    """DROP TABLE IF EXISTS _tmp_v%(visitId)d_DIASourceToObjectMatches""",
    """DROP TABLE IF EXISTS _tmp_v%(visitId)d_PredToDIASourceMatches""",
    """DROP TABLE IF EXISTS _tmp_v%(visitId)d_DIASourceToNewObject""",
]

sqlStatements = { 'mysql' : mysqlStatements }
