import iris

def select_field( workDir, pattern, fieldOptions=None ):
    """
    Read selected fields from NAME file
      * workDir is the directory in which the files are located
      * pattern is a text string common to all files to be loaded (similar)
        to a UNIX pattern. e.g. Fields_grid1_*
      * fieldOptions is a dictionary containing information about columns to load
        e.g. {'Species' : 'LPAR', 'Quantity' : 'Deposition'}
        inclusion of fieldOptions is optional
    """

    att_spec = ( {} if fieldOptions==None else fieldOptions )
    attConstraint = iris.AttributeConstraint(**att_spec)
    cubes=iris.load(workDir+'/'+pattern, attConstraint)

    return cubes
