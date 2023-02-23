import pymysql
import pandas as pd
import numpy as np
from ddb import post_process, connect_mysql, create_schema, insert_data, drop_database, drop_table, search_db, drop_row


####
# Test Dragon Descriptors from combi_lib (molecules w/ MW given, no DES)
data = pd.read_csv('combi_lib_dragon_descriptors.csv')
cols = list(data.columns)
test_properties2 = []
unit_list2 = []
molec_id_type2 = ['smiles']
mw_meta2 = False
des_status2 = False
prop_meta2 = False
# code below necessary because user must specify properties
# TODO: make insert function scan for properties by itself
for col in cols:
    if col.lower() == 'amw':
        test_properties2.append(col)
        unit_list2.append(col + ' units') # dummy units bc using dragon descriptors
    elif col.lower() != 'smiles' and col.lower() != 'mw':
        test_properties2.append(col)
        unit_list2.append(col+' units') # dummy units bc using dragon descriptors

cur, all_dbs, con = connect_mysql(host='127.0.0.1',user='root',pw='root')
drop_database(db='aab_chembddb')
create_schema(host=-1,user='root',pw='root',db='aab_chembddb')
insert_data('combi_lib_dragon_descriptors.csv',test_properties2,unit_list2,molec_id_type2,mw_meta2,des_status2,prop_meta2,'aab_chembddb')
from_form = {

    'ZM1V_id': [35],
    'ZM1V_from_val': [45],
    'ZM1V_to_val': [100.0],

    'MW': None,
    'MW_from_val': [50],
    'MW_to_val': [200],

    # 'smiles_search': None,
    # 'SMILES': ['OCC(C(=O)O)(O)C'],
    
    'molecule': None
}

# sdb = search_db('combi_lib',['search-query-molecule'],from_form)
# drop_row('combi_lib',['search-query-molecule','all'],from_form)