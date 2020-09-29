# Eval-cache data model.  This is stored as a separate table in the
# sqlite db.

import sqlite3 as s3

# How should we handle this?  We're approaching 2 levels of "cache":
# immediate (class slot) and db (persistent storage).
