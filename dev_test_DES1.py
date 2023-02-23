import pymysql
import pandas as pd
import numpy as np
from dev_des_bddb import post_process, connect_mysql, create_schema, insert_data, search_db, drop_database, drop_table, drop_row


test_properties1 = ['density','electrochemical_window','viscosity','conductivity']
test_methods1 = ['dft']
test_functionals1 = ['bp86']
test_basis1 = ['Def2-TZVP']
test_ff1 = ['UFF']
unit_list1 = ['g/cc','V','cp','mS/cm']
molec_id_type1 = ['smiles']
mw_meta1 = True
des_status1 = True
prop_meta1 = True




# Test DES experimental data
cur, all_dbs, con = connect_mysql(host='127.0.0.1',user='root',pw='legoeggo')
drop_database(db='exp_data')
create_schema(host=-1,user='root',pw='legoeggo',db='exp_data')
insert_data('DESexpDat1.csv',test_properties1,unit_list1,molec_id_type1,mw_meta1,des_status1,prop_meta1,'exp_data',methods_list = test_methods1,functionals_list = test_functionals1,basis_list = test_basis1,forcefield_list = test_ff1)
cur.execute('show tables in exp_data;')
all_db_tables = cur.fetchall()
#print(all_db_tables)

from_form = {
    
    #'density_id': [1],
    #'density_from_val': [1.0005],
    #'density_to_val': [1.15],

    #'electrochemical_window_id': [2],
    #'electrochemical_window_from_val': [2.5],
    #'electrochemical_window_to_val': [3.5],

    'viscosity_id': [3],
    'viscosity_from_val': [1.0],
    'viscosity_to_val': [10000.0],

    #'conductivity_id': [4],
    #'conductivity_from_val': [3.2],
    #'conductivity_to_val':[50.0],

    # 'HBA_MW': None,
    # 'HBA_MW_from_val': [50],
    # 'HBA_MW_to_val': [200],

    # 'HBD_MW': None,
    # 'HBD_MW_from_val': [45],
    # 'HBD_MW_to_val': [100],
    #'smiles_search': None,
    #'HBA_SMILES': ['OCC[N+](C)(C)C.[Cl-]'],
    #'HBD_SMILES': ['OCCO'],
    
    'des': None
}
sdb = search_db('exp_data',['search-query-des'],from_form)
#search_db('exp_data','download_json', from_form)
