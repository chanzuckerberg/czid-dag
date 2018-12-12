import shelve

def shelve_open_r(db_file):
    try:
        result = shelve.open(db_file, 'r')
    except:
        result = shelve.open(db_file.replace(".db", ""), 'r')
    return result
