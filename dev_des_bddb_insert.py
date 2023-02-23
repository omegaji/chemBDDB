from flask import Flask, render_template, url_for, request, redirect
import pymysql
import pandas as pd
import numpy as np
import os
import sys
from chembddb import molidentfiers
from chembddb.units import insert_unit_list,fetch_unit_list, create_unit_list, unit_converter
import json
from urllib.request import urlopen
from urllib.parse import quote
from urllib.error import HTTPError
import time
from copy import deepcopy
try:
    import pybel
except:
    from openbabel import pybel
from itertools import cycle, islice,chain

from flask import send_from_directory


all_dbs = []
app = Flask(__name__)
upload_directory = os.getcwd()
app.config['UPLOAD FOLDER']=upload_directory

def post_process(sql, from_db,meta):
    """
    Helper function to post process the results from the database

    Parameters
    ----------
    sql: str
        sql query that was created by the search function
    from_db: tuple of tuples
        results from the database
    meta: str
        specifies if Molecule or DES


    Returns
    -------
    data: pandas dataframe
        dataframe containing the processed results ready to be displayed
    columns: list of str
        list of formatted and cleaned column names
    """
    abc = 0
    if meta.lower() == 'molecule':
        if 'MW' in sql and 'property' not in sql:
            data = pd.DataFrame(list(from_db), columns=['ID','SMILES','MW'])
            columns = list(data.columns)
        elif 'MW' not in sql and 'property' not in sql:
            data = pd.DataFrame(list(from_db),columns = ['ID','SMILES'])
            columns = list(data.columns)
        else:
            data = pd.DataFrame(list(from_db), columns=['Molecule_id','SMILES','MW','Property','Value'])
            data['ID_SMI']=data['Molecule_id'].astype(str)+','+data['SMILES']+','+data['MW'].astype(str)
            data['Property']=data['Property']#+'-' +data['Method']+'('+data['Functional']+'/'+data['Basis_set']+')('+data['forcefield']+')'
            data = data[data.columns[-3:]]
            data=data.pivot_table(index='ID_SMI',columns='Property',values='Value')
            data = data.reset_index()
            data[['ID','SMILES','MW']]=data['ID_SMI'].str.split(',',expand=True)
            columns=['ID','SMILES','MW']
            for i in data.columns[1:-2]:
                columns.append(i)
            columns.pop()
            data = data[columns]
            data = data.T.drop_duplicates().T
            data.MW = data.MW.astype('float')
            columns=[c.replace('(NA/NA)','') for c in columns]
            columns=[c.replace('(na/na)','') for c in columns]
            columns=[c.replace('(NA)','') for c in columns]
            columns=[c.replace('(na)','') for c in columns]
        print("Molecule post-processing successful!")
    elif meta.lower() == 'des':
        if '_MW' in sql and 'property' not in sql:
            if 'HBA_MW' in sql and 'HBD_MW' in sql:
                data = pd.DataFrame(list(from_db),columns=['DES_ID','HBA_ID','HBA_SMILES','HBA_MW','HBD_ID','HBD_SMILES','HBD_MW','ratio'])
                columns = list(data.columns)
            elif 'HBA_MW' in sql and 'HBD_MW' not in sql:
                data = pd.DataFrame(list(from_db),columns=['DES_ID','HBA_ID','HBA_SMILES','HBA_MW','HBD_ID','HBD_SMILES','ratio'])
                columns = list(data.columns)
            elif 'HBD_MW' in sql and 'HBA_MW' not in sql:
                data = pd.DataFrame(list(from_db),columns=['DES_ID','HBA_ID','HBA_SMILES','HBD_ID','HBD_SMILES','HBD_MW','ratio'])
                columns = list(data.columns)
        elif '_MW' not in sql and 'property' not in sql:
            data = pd.DataFrame(list(from_db),columns = ['DES_ID','HBA_ID','HBA_SMILES','HBD_ID','HBD_SMILES','ratio'])
            columns = list(data.columns)
        else:
            data = pd.DataFrame(list(from_db), columns=['DES_ID','HBA_ID','HBA_SMILES','HBA_MW','HBD_ID','HBD_SMILES','HBD_MW','ratio','Property','Value'])
            data['HBA_ID_SMI'] = data['HBA_ID'].astype(str)+','+data['HBA_SMILES']
            data['HBD_ID_SMI'] = data['HBD_ID'].astype(str)+','+data['HBD_SMILES']
            data['index'] = list(range(1, len(data) + 1))
            data['ID_DES_ID_SMI'] = data['index'].astype(str)+','+data['DES_ID'].astype(str)+','+data['HBA_ID_SMI']+','+data['HBD_ID_SMI']+','+data['ratio'].astype(str)
            data['Property']=data['Property'] # TODO: add method, functionals, basis_set, forcefields to this! (refer to molecule type in chembddb)
            data = data[data.columns[-6:]]
            data=data.pivot_table(index='ID_DES_ID_SMI',columns='Property',values='Value')
            data = data.reset_index()
            data[['ID','DES_ID','HBA_ID','HBA_SMILES','HBD_ID','HBD_SMILES','ratio']] = data['ID_DES_ID_SMI'].str.split(',',expand=True)
            columns=['DES_ID','HBA_SMILES','HBD_SMILES','ratio']
            for i in data.columns[1:-4]:
                columns.append(i)
            for i in columns:
                if i == 'HBA_ID' or i == 'HBD_ID':
                    columns.remove(i)
                if i == 'ID':
                    columns.remove(i)
            columns.pop()
            data = data[columns]
            data = data.T.drop_duplicates().T
            columns=[c.replace('(NA/NA)','') for c in columns]
            columns=[c.replace('(na/na)','') for c in columns]
            columns=[c.replace('(NA)','') for c in columns]
            columns=[c.replace('(na)','') for c in columns]
        print("DES post-processing successful!")
    return data, columns

@app.route('/')
def begin():
    return redirect(url_for('connect'))

def connect_mysql(host,user,pw):
    """
    Connects user to MySQL server

    Parameters
    ----------
    host: str
        host, usually '127.0.0.1'
    user: str
        user, root
    pw: str
        password, this will depend on the user's root password


    Returns
    -------
    cur: MySQL cursor
        cursor used to query mySQL server
    all_dbs: list of str
        list of database names
    con: MySQL connection
        connection with MySQL
    
    OR

    'invalid' 'credentials' '!': str
        indicates that the user had invalid credentials
    """
    global cur, con
    try:
        con = pymysql.connect(host = host, user=user, password = pw)
        cur = con.cursor()
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        all_dbs = []
        for i in all_dbs_tup:
            if '_chembddb' in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        print('Connection successful!')
        return cur,all_dbs, con
    except:
        print('Connection failed!')
        return 'invalid','credentials', '!'

@app.route('/connect',methods=['GET','POST'])
def connect():
    """ establishes MySQL connection based on the credentials provided during setup

    Parameters
    ----------
    Returns
    -------
    cur: cursor object
        pointer to the MySQL database; used to execute all SQL queries
    """
    global cur, all_dbs, unit_list, con
    if request.method=='POST':
        cred=request.form
        cred = cred.to_dict(flat=False)
        cur,all_dbs,con = connect_mysql(host = cred['host'][0], user=cred['username'][0],pw=cred['password'][0])
        if cur == 'invalid' and all_dbs == 'credentials':
            return render_template('connect.html',err_msg='Invalid Credentials. Did not connect to MySQL.')
        else:
            print(all_dbs)
            if any('unit_list' in i for i in all_dbs):
                all_dbs.pop(all_dbs.index(('unit_list',)))
                unit_list = fetch_unit_list(cur)
            else:
                unit_list = create_unit_list(cur,con)
            return render_template('connect.html',success_msg='Connection Established',host=cred['host'][0],user=cred['username'][0],password=cred['password'][0],all_dbs=all_dbs)
    else:
        return render_template('connect.html')

@app.route('/setup',methods=['GET','POST'])
def create_schema(host=-1,user='',pw='',db=''):
    """
    Calls connect_mysql function to connect user to MySQL server
    Creates database using ChemBDDB schema

    Parameters
    ----------
    host: str default=''
        the hostname is the domain name or server name
    user: str default=''
        the username for MySQL
    pw: str default=''
        the password for MySQL
    db: str default=''
        the name of the database that needs to be set up


    Returns
    -------
    all_dbs: list of str
        list of all databases
    """
    if host != -1:
        # for python module
        b, a, c = connect_mysql(host=host,user=user,pw=pw)
        if b == 'invalid' and a == 'credentials':
            print('Invalid credentials!')
            return 'invalid credentials'
        else:
            db = db +'_chembddb'
    elif request.method=='POST':
        # for UI
        db_details=request.form
        db_details=db_details.to_dict(flat=False)
        db=db_details['dbname'][0]+'_chembddb'
    else:
        # Default landing page for setup
        all_dbs=[]
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        
        return render_template('setup.html',all_dbs=all_dbs)
    
    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result ==0:
        cur.execute('CREATE DATABASE `%s`;'%db)
        cur.execute('USE %s;'%db)
        cur.execute('CREATE TABLE `%s`.`Property`(`id` INT NOT NULL AUTO_INCREMENT,`Property_str` VARCHAR(100) NOT NULL UNIQUE,`Unit` VARCHAR(100) NOT NULL,PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Molecule` (`id` INT NOT NULL AUTO_INCREMENT,`SMILES` VARCHAR(300) DEFAULT \'NONE\', `Standard_InChI` VARCHAR(400) DEFAULT \'NONE\',`Standard_InChI_Key` VARCHAR(100) DEFAULT \'NONE\',`CAS_Registry_Number` VARCHAR(200) DEFAULT \'NONE\',`IUPAC_Name` VARCHAR(400) DEFAULT \'NONE\',`Other_name` VARCHAR(1000) DEFAULT \'NONE\',`Chemical_Formula` VARCHAR(100) DEFAULT \'NONE\',`MW` FLOAT NOT NULL, PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Model`(`id` INT NOT NULL AUTO_INCREMENT,`method_name` VARCHAR(100) NOT NULL UNIQUE,`options` VARCHAR(500),PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Functional`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Basis_set`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Forcefield`(`id` INT NOT NULL AUTO_INCREMENT,`name` VARCHAR(100) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        # DES entries cannot have identical hba and hbd ids from molecule table
        cur.execute('CREATE TABLE `%s`.`DES` (`id` INT NOT NULL AUTO_INCREMENT,`HBA_id` INT,`HBD_id` INT,`ratio` FLOAT, PRIMARY KEY(`id`),FOREIGN KEY(`HBA_id`) REFERENCES `Molecule`(`id`) ON DELETE CASCADE,FOREIGN KEY(`HBD_id`) REFERENCES `Molecule`(`id`) ON DELETE CASCADE,CHECK(`HBA_id`<>`HBD_id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Value` (`id` INT NOT NULL AUTO_INCREMENT,`num_value` FLOAT NOT NULL,`model_id` INT,`Property_id` INT NOT NULL, `functional_id` INT, `basis_id` INT, `forcefield_id` INT, `DES_id` INT, `molecule_id` INT, PRIMARY KEY(`id`));'%db)
        cur.execute('CREATE TABLE `%s`.`Configuration`(`id` INT DEFAULT 0,`conf` VARCHAR(200) DEFAULT \'NONE\',`unit_dict` VARCHAR(500) DEFAULT \'NONE\',PRIMARY KEY (`id`));'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk0` FOREIGN KEY (`model_id`) REFERENCES `Model`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk1` FOREIGN KEY (`property_id`) REFERENCES `Property`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk2` FOREIGN KEY (`functional_id`) REFERENCES `Functional`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk3` FOREIGN KEY (`basis_id`) REFERENCES `Basis_set`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk4` FOREIGN KEY (`forcefield_id`) REFERENCES `Forcefield`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk5` FOREIGN KEY (`DES_id`) REFERENCES `DES`(`id`) on DELETE CASCADE;'%db)
        cur.execute('ALTER TABLE `%s`.`Value` ADD CONSTRAINT `Value_fk6` FOREIGN KEY (`molecule_id`) REFERENCES `Molecule`(`id`) on DELETE CASCADE;'%db)
        cur.execute('SHOW TABLES FROM `%s`;'%db)
        all_dbs_tup=cur.fetchall()
        all_dbs=[]
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        print('The database has been created!')
        if host == -1:
            # successful creation for UI
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,success_msg='The database has been created.')
        else:
            # successful creation for python module
            return 'Success'
    else:
        if host == -1:
            # error handling for UI
            return render_template('setup.html',dbname=db,all_dbs=all_dbs,err_msg='Database already exists.')
        else:
            # error handling for python module
            return 'Failed! Database already exists.'

@app.route('/temp_insert',methods=['GET','POST'])
def temp_insert():
    """
    Inserts csv data into database

    """
    global all_dbs, db, data_file, cur, data,mol_ids,con,snapshot,prop_type,prop_store,sim_status,mw_cols,mw_meta
    global des_status,methods_list,functionals_list,basis_list,forcefield_list, molecule_identifiers, molecule_identifiers_cols
    global pybel_identifiers, molecule_id, mw_cols
    mi_cols = []
    cur.execute('show databases;')
    all_dbs_tup = cur.fetchall()
    all_dbs = []
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    if request.method == 'POST' and 'upload_data' in request.form:
        config_options = request.form
        config_options=config_options.to_dict(flat=False)
        db = config_options['dbname'][0]
        files = request.files
        cur.execute('USE {}_chembddb;'.format(db))
        cur.execute('SELECT ID, conf from Configuration')
        confs = cur.fetchall()
        conf = False
        if len(confs) > 0:
            conf = True
        present = []
        identifiers = ['SMILES','Standard_Inchi_Key','Standard_Inchi','Chemical_Formula','IUPAC_Name','Other_name','CAS_Registry_Number']
        snapshot = []
        all_dbs_tup=cur.fetchall()
        all_dbs=[]
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        data_file = files['data_file']
        if data_file is not None:
            data_file.seek(0)
            data = pd.read_csv(data_file)
        
        cols = list(data.columns)
        pybel_identifiers = {'smiles':'smi','standard_inchi_key':'inchikey','standard_inchi':'inchi'}
        #prop_val = []
        #unit_val = []
        molec_id = []
        id_col_name = []
        mol_ids = {}
        for id in identifiers:
            cur.execute('SELECT COUNT(DISTINCT '+ id +') FROM MOLECULE;')
            counts = cur.fetchall()
            if counts[0][0] > 1:
                present.append(id)
        if present !=[]:
            snapshot = [', '.join(p for p in present)]
            snapshot = [snapshot[0].split(',')]
            # check which properties are present in the database

            cur.execute('SELECT Property_str, Unit from Property;')
            properties = cur.fetchall()
            properties = ', '.join(i[0]+'('+i[1]+')' for i in properties)
            properties.replace(', na','')
            properties.replace('na,','')
            properties.replace('na','')

            properties = properties.split(',')
            snapshot.append(properties)
            
            cur.execute('SELECT method_name from Model;')
            methods = cur.fetchall()
            methods = ', '.join(i[0] for i in methods)
            methods = methods.replace(', na','')
            methods = methods.replace ('na,','')
            methods = methods.replace('na','')

            methods = methods.split(',')
            snapshot.append(methods)

            cur.execute('select name FROM Functional;')
            functionals = cur.fetchall()
            functionals = ', '.join(i[0] for i in functionals)
            functionals = functionals.replace(', na','')
            functionals = functionals.replace('na,','')
            functionals = functionals.replace('na','')

            functionals = functionals.split(',')
            snapshot.append(functionals)

            cur.execute('SELECT name FROM Basis_set;')
            basis = cur.fetchall()
            basis = ', '.join(i[0] for i in basis)
            basis = basis.replace(', na','')
            basis = basis.replace('na,','')
            basis = basis.replace('na','')

            basis = basis.split(',')
            snapshot.append(basis)

            cur.execute('SELECT name FROM Forcefield')
            forcefield = cur.fetchall()
            forcefield = ', '.join(i[0] for i in forcefield)
            forcefield = forcefield.replace(', na','')
            forcefield = forcefield.replace('na,','')
            forcefield = forcefield.replace('na','')
            forcefield = forcefield.replace(' ','')
            forcefield = forcefield.split(',')
            snapshot.append(forcefield)
   
        else:
            snapshot = False
        
        
        if data_file.filename.rsplit('.',1)[1]!='csv':
            db.replace('_chembddb','')
            db = db.replace('_',' ')
            return render_template('temp_insert.html',all_dbs=all_dbs,title=db,err_msg='No data file provided or incorrect file format. (csv requred)')
        else:
            des_status = request.form['des_status']
            prop_store = request.form['prop_store']
            #print(prop_store)
            cols = []
            for i in data.columns:
                cols.append(i.replace(' ','_'))
            print(data.columns)
            #print(i)
            data.columns = cols
            if prop_store == 'single_col':
                prop_store = True
            else:
                prop_store = False
            sim_status = request.form['sim_status']
            if sim_status == 'yes_sim':
                sim_status = True
                methods_list = True
                functionals_list = True
                basis_list = True
                forcefield_list = True
            else:
                sim_status = False
                methods_list = False
                functionals_list = False
                basis_list = False
                forcefield_list = False
            
            mw_meta = request.form['mw_meta']
            if mw_meta == 'yes_mw':
                mw_meta = True
            else:
                mw_meta = False

            if des_status == 'molecule':
                des_status=False
                return render_template('temp_insert.html',all_dbs=all_dbs,data_validated=True,des_status=False,sim_status=sim_status,prop_store=prop_store,mw_meta=mw_meta,cols = list(data.columns),conf=conf,snapshot=snapshot,snapshot_cols = ['Molecule Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields'])
            else:
                des_status=True
                return render_template('temp_insert.html',all_dbs=all_dbs,data_validated=True,des_status=True,sim_status=sim_status,prop_store=prop_store,mw_meta=mw_meta,cols=list(data.columns),conf=conf,snapshot=snapshot,snapshot_cols = ['HBA_Identifiers','HBD_Identifiers','Porperties (Units)','Methods','Functionals','Basis Sets','Forcefields'])
    elif request.method == 'POST' and ('config' in request.form or 'use-config' in request.form):
        
        config_options = request.form
        config_options = config_options.to_dict(flat=False)
        molecule_identifiers_cols = []
        molecule_identifiers = []
        mw_cols = []
        cols = list(data.columns)
        
        for key,value in config_options.items():
            if 'col' in key:
                molecule_identifiers_cols.append(value)
            elif '_id' in key:
                if value != ['cols_0A','cols_0B'] and value != ['MW']:
                    molecule_identifiers.append(value)
            elif 'mw' in key:
                if value != ['cols_0A_MW','cols_0B_MW']:
                    mw_cols.append(value)
        #molecule_identifiers = list(chain(*molecule_identifiers))
        

        remaining_cols = []
        for i in cols:
            if i not in molecule_identifiers:
                remaining_cols.append(i)
        if des_status == True:
            snapshot_cols = ['HBA_Identifiers','HBD_Identifiers','ratio','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields']
        else:
            snapshot_cols = ['Molecule_Identifiers','Properties (Units)','Methods','Functionals','Basis Sets','Forcefields']
        
        
        return render_template('temp_insert.html',props=remaining_cols,prop_length=len(remaining_cols),all_dbs=all_dbs,title=db,snapshot=snapshot,snapshot_cols=snapshot_cols,prop_store=prop_store,des_status=des_status,sim_status=sim_status,mw_meta=mw_meta)

    elif request.method == 'POST' and ('meta-data' in request.form or 'download-submit' in request.form):
        print("IN THE FORM")
        print(request.form)
        meta_data = request.form
        meta_data = meta_data.to_dict(flat=False)
        
        property_list = meta_data['2_prop']
        unit_list = meta_data['2_unit']
        if sim_status == True:
            simu_data_columns = meta_data['sim_id_0']
        if prop_store == True:
            property_columns = meta_data['prop_store_id_0']
        else:
            property_columns = meta_data['prop_id_0']
        mol_ids={}

        if des_status == True:
            full_id_col = data[molecule_identifiers[0][0]].values.tolist() + data[molecule_identifiers[0][1]].values.tolist()
        else:
            full_id_col = data[molecule_identifiers[0]].values.tolist()
        get_id = full_id_col
        molec_id = []
        for id in get_id:
            if id not in molec_id:
                molec_id.append(id)
            else:
                pass
        if des_status == True:
            mol_ids.update({molecule_identifiers_cols[0][0]:molec_id})
        else:
            mol_ids.update({molecule_identifiers_cols[0][0]:molec_id})
        # Add Properties to Property Table
        cur.execute('SELECT Property_str FROM Property;')
        old_properties = cur.fetchall()
        old_properties = tuple(tuple(x) for x in old_properties)
        cur.execute("SELECT unit FROM Property;")
        old_units = cur.fetchall()
        old_units = tuple(tuple(x) for x in old_units)
        new_properties = []
        new_units = []
        for i in property_list:
            if i not in old_properties:
                new_properties.append(i)
        new_properties = tuple(new_properties)
        for i in unit_list:
            if i not in old_units:
                new_units.append(i)
        new_units = tuple(new_units)
        for i in range(len(new_properties)):
            cur.execute('INSERT INTO Property(Property_str,unit) VALUE ("%s","%s");'%(new_properties[i],new_units[i]))
        # now for methods, functionals, basis sets, forcefields
        if methods_list != False:
            methods_list = data[simu_data_columns[0]]
            cur.execute('SELECT method_name from Model')
            old_methods = cur.fetchall()
            old_methods = tuple(tuple(x) for x in old_methods)
            new_methods = []
            for i in methods_list:
                if i not in old_methods:
                    if i not in new_methods:
                        new_methods.append(i)
            new_methods = tuple(new_methods)
            for i in new_methods:
                cur.execute('INSERT INTO Model(method_name) VALUE("%s")'%i)
            cur.execute('SELECT id, method_name from Model')
            method_info = cur.fetchall()
            method_info = tuple(tuple(x) for x in method_info)
            method_id = []
            methods = {}
            for x in method_info:
                methods.update({x[0]: x[1]})
            for k, v in methods.items():
                name_col = simu_data_columns[0]
                for i in data[name_col]:
                    if i == v:
                        method_id.append(k)

        if functionals_list != False:
            functionals_list = data[simu_data_columns[1]]
            cur.execute('SELECT name from Functional;')
            old_functionals = cur.fetchall()
            old_functionals = tuple(tuple(x) for x in old_functionals)
            new_functionals = []
            for i in functionals_list:
                if i not in old_functionals:
                    if i not in new_functionals:
                        new_functionals.append(i)
            new_functionals = tuple(new_functionals)
            for i in new_functionals:
                cur.execute('INSERT INTO Functional(name) VALUE("%s")'%i)
            cur.execute('SELECT id, name from Functional')
            functional_info = cur.fetchall()
            functional_info = tuple(tuple(x) for x in functional_info)
            functional_id = []
            functionals = {}
            for x in functional_info:
                functionals.update({x[0]: x[1]})
            for k, v in functionals.items():
                name_col = simu_data_columns[1]
                for i in data[name_col]:
                    if i == v:
                        functional_id.append(k)
        
        if basis_list != False:
            basis_list = data[simu_data_columns[2]]
            cur.execute('SELECT name from Basis_set;')
            old_basis = cur.fetchall()
            old_basis = tuple(tuple(x) for x in old_basis)
            new_basis = []
            for i in basis_list:
                if i not in old_basis:
                    if i not in new_basis:
                        new_basis.append(i)
            new_basis = tuple(new_basis)
            for i in new_basis:
                cur.execute('INSERT INTO Basis_set(name) VALUE("%s")'%i)
            cur.execute('SELECT id, name from Basis_set')
            basis_info = cur.fetchall()
            basis_info = tuple(tuple(x) for x in basis_info)
            basis_id = []
            basis = {}
            for x in basis_info:
                basis.update({x[0]: x[1]})
            for k, v in basis.items():
                name_col = simu_data_columns[2]
                for i in data['basis']:
                    if i == v:
                        basis_id.append(k)
        
        if forcefield_list != False:
            forcefield_list = data[simu_data_columns[3]]
            cur.execute('SELECT name from Forcefield;')
            old_forcefield = cur.fetchall()
            old_forcefield = tuple(tuple(x) for x in old_forcefield)
            new_forcefield = []
            for i in forcefield_list:
                if i not in old_forcefield:
                    if i not in new_forcefield:
                        new_forcefield.append(i)
            new_forcefield = tuple(new_forcefield)
            for i in new_forcefield:
                cur.execute('INSERT INTO Forcefield(name) VALUE("%s")'%i)
            cur.execute('SELECT id, name from Forcefield')
            forcefield_info = cur.fetchall()
            forcefield_info = tuple(tuple(x) for x in forcefield_info)
            forcefield_id = []
            forcefield = {}
            for x in forcefield_info:
                forcefield.update({x[0]: x[1]})
            for k, v in forcefield.items():
                name_col = simu_data_columns[3]
                for i in data['forcefield']:
                    if i == v:
                        forcefield_id.append(k)
        return render_template('temp_insert.html',all_dbs=all_dbs,init=True,success_msg='Success')
    
    elif request.method == 'POST' and ('final_md' in request.form):
        # Add Molecules to Molecule Table
        new_entries = []
        row = []
        if mw_meta == False:
            mw_flag = False
            for k in mol_ids.keys():
                for val in mol_ids[k]:
                    if k.lower() in pybel_identifiers.keys():
                        try:
                            m = pybel.readstring(pybel_identifiers[k.lower()],val)
                            if k.lower() == 'smiles':
                                iden = m.write('can').strip()
                                if mw_flag == False:
                                    mw = m.molwt
                                    mw = round(mw,3)
                                    mw_flag = True
                            else:
                                iden = m.write(pybel_identifiers[k.lower()]).strip()
                                if mw_flag == False:
                                    mw = m.molwt
                                    mw = round(mw,3)
                                    mw_flag = True
                            row.append((iden,mw))
                            del(mw)
                            del(iden)
                            mw_flag = False
                        except:
                            for mol in range(len(data)):
                                db = db.replace('_',' ')
                                print('Invalid Smiles on row number '+ str(mol) + '!')
                    else:
                            row.append(val)
                            if mw_flag == False:
                                mw_flag = True
                                url = 'http://cactus.nci.nih.gov/chemical/structure/'
                                try:
                                    url = url + val + '/mw'
                                    ans = urlopen(url).read().decode('utf8')
                                except HTTPError:
                                    mw = NULL
            new_entries.append(row)
            new_entries = tuple(tuple(x) for x in new_entries)
        else:
            # if user used library generator that gives MW values
            mw_data = {}
            if des_status == True:
              mw = []
              mw = data[mw_cols[0]].values.tolist() + data[mw_cols[1]].values.tolist()
            else:
                mw = data[mw_cols[0]].values.tolist()
            print(mw)
            for i in range(len(get_id)):
                mw_data[get_id[i]] = mw[i]
            for i in mw_data.keys():
                new_entries.append((i,mw_data[i]))
        if mw_meta == True:
            # insert information into Molecule Table
            cur.execute("SELECT "+ ''.join(i.lower()+',' for i in mol_ids.keys())+"MW from Molecule")
            molecules = cur.fetchall()
            molecules = tuple(tuple(x) for x in molecules)
            required_entries = list(set(new_entries)-set(molecules))
            if len(mol_ids.keys())>1:
                mol_q = ''.join(i+',' for i in mol_ids.keys())
                vals = ''.join('%s,' for i in range(len(mol_ids)+1))
            else:
                mol_q = list(mol_ids.keys())[0] + ','
                vals = '%s,%s,'
            for x in required_entries:
                can_smiles = pybel.readstring('smi',x[0]).write('can').strip()
                print(can_smiles)
                cur.execute('INSERT INTO Molecule('+mol_q+'MW) VALUE('+vals[:-1]+')',(can_smiles,x[1]))
            print('Molecules done!')
        else:
            cur.execute("SELECT "+ ''.join(i.lower()+',' for i in mol_ids.keys())+"MW from Molecule")
            molecules = cur.fetchall()
            molecules = tuple(tuple(x) for x in molecules)
            required_entries = list(set(new_entries)-set(molecules))
            
            if len(mol_ids.keys())>1:
                mol_q = ''.join(i+',' for i in mol_ids.keys())
                vals = ''.join('%s,' for i in range(len(mol_ids)+1))
            else:
                mol_q = list(mol_ids.keys())[0] + ','
                vals = '%s,%s,'

            for x in required_entries[0]:
                can_smiles = pybel.readstring('smi',x[0]).write('can').strip()
                cur.execute('INSERT INTO Molecule('+mol_q+'MW) VALUE('+vals[:-1]+')',(can_smiles,x[1]))
            print('Molecules done!')

        cur.execute('SELECT id,Property_str from Property')
        all_props = cur.fetchall()
        prop_id = dict(map(reversed,all_props)) # reversed so that keys are the names of the properties and values are the id numbers
        mol_q = 'ID,'+ list(mol_ids.keys())[0]
        cur.execute("SELECT "+mol_q+" from Molecule")
        all_mols = cur.fetchall()
        molecule_id = dict(map(reversed,all_mols))
        insert_df = pd.DataFrame()
        if prop_store == True:
            prop_val = data[property_columns[0]].values.tolist()
            prop_type = data[property_columns[1]].values.tolist()
            property_id = []
            for i in prop_type:
                property_id.append(prop_id[i]) # list that will be put into insert_df. Uses the name of the property as key to get id as value.
            insert_df['property_id'] = property_id
            insert_df['value'] = prop_val
        else:
            descriptors = data[property_list].melt() # each molecule comes up 1 time for each property
            prop_type = descriptors['variable'].values.tolist()
            prop_val = descriptors['value'].values.tolist()
            property_id = []
            for i in prop_type:
                property_id.append(prop_id[i]) # list that will be put into insert_df. Uses the name of the property as key to get id as value.
            insert_df['property_id'] = property_id
            insert_df['value'] = prop_val
        if des_status == False:
            temp_df = pd.DataFrame({'molecule':get_id}) 
            temp_df['mol_ids'] = temp_df['molecule'].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
            mol_ids_list = temp_df['mol_ids'].values.tolist()
            insert_df['molecule_id'] = list(islice(cycle(mol_ids_list),len(insert_df)))
        else:

            des = {}
            distinct_des_dict = {}
            #print(data.head())
            des_col_info = np.concatenate(molecule_identifiers).flat
            #print(des_col_info)
            #print(data.loc[:,des_col_info].head())

            hba_list = data[des_col_info[0]].values.tolist()
            hbd_list = data[des_col_info[1]].values.tolist()
            ratio_list = data[des_col_info[2]].values.tolist()
            for i in range(len(data)):
                des.update({i:tuple([hba_list[i],hbd_list[i],ratio_list[i]])})

            distinct_des = list(set(des.values()))

            for i in range(len(distinct_des)):
                distinct_des_dict.update({i:distinct_des[i]})
            
            distinct_des_df = pd.DataFrame(distinct_des_dict)
            distinct_des_df = distinct_des_df.T
            
            distinct_des_df[0] = distinct_des_df[0].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
            distinct_des_df[1] = distinct_des_df[1].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
            
            insert_df['HBA'] = hba_list
            #print(insert_df)
            insert_df['HBA_id']=insert_df['HBA'].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
            insert_df['HBD'] = hbd_list
            insert_df['HBD_id']=insert_df['HBD'].apply(lambda a: molecule_id[pybel.readstring('smi',a).write('can').strip()])
            insert_df = insert_df.drop(['HBA','HBD'],axis=1)
            insert_df['ratio'] = ratio_list
            

            for i in range(len(distinct_des_df)):
                hba = distinct_des_df[0][i]
                hbd = distinct_des_df[1][i]
                ratio = distinct_des_df[2][i]
                cur.execute('INSERT INTO DES(HBA_id,HBD_id,ratio) VALUES ("%s","%s",%s);'%(hba,hbd,ratio))

            cur.execute('SELECT id, HBA_id, HBD_id, ratio from DES')
            all_des = cur.fetchall()
            des_key = []
            des_val = []
            for i in all_des:
                des_key.append(i[0])
                des_val.append(str(i[1])+str(i[2])+str(i[3]))
            des_dict = dict(map(reversed,zip(des_key,des_val)))
            des_temp = []
            for i in range(len(insert_df)):
                des_entry = str(insert_df['HBA_id'][i])+str(insert_df['HBD_id'][i])+str(round(insert_df['ratio'][i],6))
                des_temp.append(des_entry)
            insert_df['des_temp'] = des_temp
            insert_df['DES_id'] = insert_df['des_temp'].apply(lambda a: des_dict[a])
            insert_df = insert_df.drop(['HBA_id','HBD_id','ratio','des_temp'],axis=1)
            print('DES done!')
        
        mol_or_des_only = False
        if methods_list != False:
            insert_df['model_id'] = method_id
        else:
            mol_or_des_only = True
        if functionals_list != False:
            insert_df['functional_id'] = functional_id
        else:
            mol_or_des_only = True
        if basis_list != False:
            insert_df['basis_id'] = basis_id
        else:
            mol_or_des_only = True
        if forcefield_list != False:
            insert_df['forcefield_id'] = forcefield_id
        else:
            mol_or_des_only = True
        
        # Now, taking data to insert into values table
        print('What is being entered into Values Table:')
        print(insert_df)
        if mol_or_des_only == True:
            print('Data for Molecule/DES is being added!')
            if des_status == False:
                mol_cols = ['property_id','value','molecule_id']
                insert_df = insert_df[mol_cols]
                insert_tuples = insert_df.to_records(index=False)
                insert_tuples = insert_tuples.values.tolist() # execute many needs sequence of sequences to work
                print('Inserting data into Values table...')
                start = time.time()
                cur.executemany('INSERT INTO Value(property_id,num_value,molecule_id) VALUES(%s,%s,%s);',(insert_tuples))
                stop = time.time()
                print('This took '+str(round((stop-start),4))+ ' seconds!')
            else:
                des_cols = ['property_id','value','DES_id']
                insert_df = insert_df[des_cols]
                insert_tuples = insert_df.to_records(index=False)
                insert_tuples = insert_tuples.values.tolist() # execute many needs sequence of sequences to work
                print('Inserting data into Values table...')
                start = time.time()
                cur.executemany('INSERT INTO Value(property_id,num_value,des_id) VALUES(%s,%s,%s);',(insert_tuples))
                stop = time.time()
                print('This took '+str(round((stop-start),4))+ ' seconds!')
        else:
            print('Data for Models, Functionals, Basis Sets, and Forcefields is being added!')
            if des_status == False:
                insert_tuples = insert_df.to_records(index=False)
                insert_tuples = insert_tuples.values.tolist() # execute many needs sequence of sequences to work
                print('Inserting data into Values table...')
                start = time.time()
                cur.executemany('INSERT INTO Value(property_id,num_value,molecule_id,model_id,functional_id,basis_id,forcefield_id) VALUES(%s,%s,%s,%s,%s,%s,%s);',(insert_tuples))
                stop = time.time()
                print('This took '+str(round((stop-start),4))+ ' seconds!')
            else:
                insert_tuples = insert_df.to_records(index=False)
                insert_tuples = insert_tuples.values.tolist() # execute many needs sequence of sequences to work
                print('Inserting data into Values table...')
                start = time.time()
                cur.executemany('INSERT INTO Value(property_id,num_value,des_id,model_id,functional_id,basis_id,forcefield_id) VALUES(%s,%s,%s,%s,%s,%s,%s);',(insert_tuples))
                stop = time.time()
                print('This took '+str(round((stop-start),4))+ ' seconds!')
        con.commit()
        print('Values inserted successfully!')
        return render_template('temp_insert.html',val_submit = True,all_dbs=all_dbs,msg="Values submitted successfully!")
    else:
        # default landing page
        return render_template('temp_insert.html',all_dbs=all_dbs,init='True',snapshot='')
        
    

    
