#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file   ostap/utils/monkey_patch.py  
# Identify all ``monkey-patch'' occurences for the given module/file 
# @author Vanya Belyaev Ivan.Belyaev@itep.ru
# @date   2016-12-02
# =============================================================================
""" Identify all ``monkey-patch'' occurences for the given module/file 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'monkey_patch' , ## all ``monkey-patch'' occurenced for the given module/file 
    )
# =============================================================================
import ast, inspect, os, sys, zipfile
# =============================================================================
## Scans a Python file via AST and finds attribute assignments to external classes/modules 
# (i.e., potential monkey patching), excluding classes defined locally within the file.
# 
# :param target: Path to file (str), a module object, or None. If None, scans the file of the caller.
# 
def monkey_patch ( target = None , logger = None ):
    """ Scans a Python file via AST and finds attribute assignments to external classes/modules 
    (i.e., potential monkey patching), excluding classes defined locally within the file.
    
    :param target: Path to file (str), a module object, or None. 
                   If None, scans the file of the caller.
    """
    
    file_path = None
    zip_archive_path = None
    zip_inner_path = None
    
    # 1. Determine the file path based on the provided target argument
    if target is not None:
        if isinstance ( target , str ):
            if target.endswith ( '.zip' ) and os.path.exists ( target ):
                print(f"⚠️ Error: '{target}' is a ZIP archive, please specify a .py file inside it or use a module.")
                return
            if '.zip/' in target or '.egg/' in target:
                parts = target.replace ( '.egg/' , '.zip/' ).split ( '.zip/' )
                zip_archive_path = parts[0] + '.zip'
                zip_inner_path = parts[1]
                if not os.path.exists ( zip_archive_path ):
                    print(f"⚠️ Error: ZIP archive not found at '{zip_archive_path}'.")
                    return
            else:
                if not target.endswith ( '.py' ) or not os.path.exists ( target ) or not os.path.isfile ( target ):
                    print(f"⚠️ Error: '{target}' is not a valid Python (.py) file path or does not exist.")
                    return
                file_path = target
        elif inspect.ismodule ( target ):
            file_path = getattr ( target ,  '__file__'   , None )
        else:
            mod_name  = getattr ( target ,  '__module__' , None )
            if mod_name and mod_name in sys.modules:
                file_path = getattr ( sys.modules [mod_name] , '__file__' , None )
    else:
        # Fall back to the caller's file if no target is specified
        caller_frame = inspect.currentframe().f_back
        file_path    = caller_frame.f_globals.get ( '__file__' , None )

    if file_path and ( '.zip/' in file_path or '.egg/' in file_path ):
        parts = file_path.replace ( '.egg/' , '.zip/' ).split ( '.zip/' )
        zip_archive_path = parts[0] + '.zip'
        zip_inner_path = parts[1]
        file_path = None

    if not zip_archive_path and ( not file_path or not os.path.exists ( file_path ) ):
        print(f"⚠️ Could not determine file path for analysis (path: {file_path}).")
        return

    # Read source code (either from file or ZIP)
    source_code = ""
    display_name = ""

    if zip_archive_path and zip_inner_path:
        display_name = f"{os.path.basename(zip_archive_path)} -> {zip_inner_path}"
        try:
            with zipfile.ZipFile ( zip_archive_path , 'r' ) as z:
                with z.open ( zip_inner_path ) as f:
                    source_code = f.read().decode ( 'utf-8' )
        except Exception as e:
            print(f"⚠️ Error reading file from ZIP archive: {e}")
            return
    else:
        file_path = os.path.abspath ( file_path )
        display_name = os.path.basename ( file_path )
        try:
            with open ( file_path , 'r' , encoding = 'utf-8' ) as f:
                source_code = f.read()
        except Exception as e:
            print(f"⚠️ Error reading file '{display_name}': {e}")
            return

    try:
        tree = ast.parse ( source_code , filename = display_name )
    except SyntaxError as e:
        print(f"⚠️ Syntax error while parsing '{display_name}': {e}")
        return

    # 2. Collect all locally defined class names in this file
    local_classes = set ()
    for node in ast.walk(tree):
        if isinstance ( node , ast.ClassDef ):
            local_classes.add ( node.name )

    patches = []

    # 3. Helper function to recursively process nodes and handle loops (e.g., for d in (...))
    def process_node(node, active_loop_vars=None):
        if active_loop_vars is None:
            active_loop_vars = {}

        if isinstance(node, ast.For):
            loop_var_name = None
            if isinstance(node.target, ast.Name):
                loop_var_name = node.target.id
            
            iterable_items = []
            if isinstance(node.iter, (ast.Tuple, ast.List)):
                iterable_items = [ast.unparse(elt) for elt in node.iter.elts]
            
            new_loop_vars = active_loop_vars.copy()
            if loop_var_name and iterable_items:
                new_loop_vars[loop_var_name] = iterable_items

            for subnode in node.body:
                process_node(subnode, new_loop_vars)
                
        elif isinstance(node, ast.Assign):
            for target_node in node.targets:
                if isinstance(target_node, ast.Attribute):
                    target_expr = ast.unparse ( target_node )
                    
                    base_node = target_node.value
                    while isinstance ( base_node, ast.Attribute ) :
                        base_node = base_node.value
                    
                    base_name = None
                    if isinstance ( base_node , ast.Name ):
                        base_name = base_node.id
                    
                    # Skip if the base class/object is defined locally in this file or is 'self'
                    if base_name in local_classes or base_name == 'self' : 
                        continue
                        
                    # If base is a loop variable, expand it for each item in the loop iterable
                    if base_name in active_loop_vars:
                        for item in active_loop_vars[base_name]:
                            expanded_expr = target_expr.replace ( f"{base_name}.", f"{item}." )
                            patches.append ( { 'line': node.lineno, 'target': expanded_expr } )
                    else:
                        patches.append ( { 'line': node.lineno, 'target': target_expr } )
        else:
            for child in ast.iter_child_nodes(node):
                process_node(child, active_loop_vars)

    # Run the recursive AST processor
    process_node ( tree )

    # 4. Print results    
    if not logger :
        from ostap.logger.logger import getLogger
        logger = getLogger ( "ostap.utils.monkey_patch" )

    if not logger :
        from ostap.logger.logger import getLogger
        logger = getLogger ( "ostap.utils.monkey_patch" )

    if not patches :
        from ostap.logger.colorized import allright
        message = f"✅ No monkey patches in '{display_name}'"        
        logger.info ( allright ( message ) ) 
    else :
        title = f"🚨 Monkey patches in '{display_name}'" 
        rows  = [ ( '#line' , 'Target' ) ]
        for p in patches:
            row = '%d' % p['line'] , p['target']
            rows.append ( row )
        import ostap.logger.table as T
        table = T.table ( rows , title = title , prefix = '# ' , alignment = 'cw' )
        logger.info ( '%s:\n%s' % ( title , table ) )
        
# =============================================================================
if '__main__' == __name__ :

    from ostap.logger.logger import getLogger
    logger = getLogger ( "ostap.utils.monkey_patch" )
    
    from ostap.utils.docme import docme 
    docme ( __name__ ) 

    for f in sys.argv : monkey_patch ( f , logger = logger )
        
# =============================================================================
##                                                                      The END 
# =============================================================================
