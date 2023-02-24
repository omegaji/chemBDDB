from flask import Flask, render_template, url_for, request, redirect,jsonify
from werkzeug.utils import secure_filename
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
from itertools import cycle, islice

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
    except Exception as e:
        print('Connection failed!')
        print(e)
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
    global cur, all_dbs, unit_list, con, cred
    if request.method=='POST':
        cred=request.form
        cred = cred.to_dict(flat=False)
        cur,all_dbs,con = connect_mysql(host = cred['host'][0], user=cred['username'][0],pw=cred['password'][0])
        print(cred)
        if cur == 'invalid' and all_dbs == 'credentials':
            
            return render_template('connect.html',err_msg='Invalid Credentials. Did not connect to MySQL.')
        else:
            if any('unit_list' in i for i in all_dbs):
                all_dbs.pop(all_dbs.index(('unit_list',)))
                unit_list = fetch_unit_list(cur)
            else:
                unit_list = create_unit_list(cur,con)
            
            return render_template('setup.html',success_msg='Connection Established',host=cred['host'][0],user=cred['username'][0],password=cred['password'][0],all_dbs=all_dbs)
    else:
        try:
            if cur!=None and cur!='invalid':
                print(cur)
                return render_template('setup.html',success_msg='Connection Established',host=cred['host'][0],user=cred['username'][0],password=cred['password'][0],all_dbs=all_dbs)

            else:
                return render_template('connect.html')
        except Exception as e:
            print("in the except part",e)
            return render_template('connect.html')

@app.route('/setup',methods=['POST'])
def create_schema(host=-1,user='',pw='',db=''):
    print(host,user,pw,db)
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
       
        db=db_details['dbname']+'_chembddb'
        print(" in post")
        print(db_details)
    else:
        # Default landing page for setup
        all_dbs=[]
        cur.execute('show databases;')
        all_dbs_tup=cur.fetchall()
        print("all dbsssss")
        print(all_dbs)
        for i in all_dbs_tup:
            if '_chembddb' in i[0] and 'unit_list' not in i[0]:
                m=i[0]
                all_dbs.append((m[:-9],))
        
        return render_template('setup.html',all_dbs=all_dbs)
    
    all_dbs=[]
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    print("fetched dbs")
    for i in all_dbs_tup:
        if '_chembddb' in i[0] and 'unit_list' not in i[0]:
            m=i[0]
            all_dbs.append((m[:-9],))
    all_dbs.append((db[:-9],))
    cur.execute('USE INFORMATION_SCHEMA')
    result=cur.execute('SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME=\'%s\''%db)
    if result ==0:
        print(" about to create the db")
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
        # all_dbs_tup=cur.fetchall()
        # all_dbs=[]
        # for i in all_dbs_tup:
        #     if '_chembddb' in i[0] and 'unit_list' not in i[0]:
        #         m=i[0]
        #         all_dbs.append((m[:-9],))
        print('The database has been created!')
        if host == -1:
            # successful creation for 
            
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
        print(request.form)
        print(meta_data)
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
        new_properties = set()
        new_units = []
        for i in property_list:
            if i not in old_properties:
                new_properties.add(i)
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
        #return render_template('temp_insert.html',all_dbs=all_dbs,init=True,success_msg='Success')
    
    # elif request.method == 'POST' and ('final_md' in request.form):
        # Add Molecules to Molecule Table
        new_entries = []
        row = []
        if mw_meta == False:
            mw_flag = False
            for k in mol_ids.keys():
                for val in mol_ids[k]:
                    if k.lower() in pybel_identifiers.keys():
                        try:
                            
                            m = pybel.readstring(pybel_identifiers[k.lower()],val[0])
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
            t=temp_df["molecule"].values
           
            temp_df['mol_ids'] = temp_df['molecule'].apply(lambda a: molecule_id[pybel.readstring('smi',a[0]).write('can').strip()])
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
                insert_tuples = insert_tuples.tolist() # execute many needs sequence of sequences to work
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
        
    
@app.route('/search',methods=['GET','POST'])
def search():
    global all_dbs
    return render_template('search.html',all_dbs=all_dbs)

@app.route('/search_db',methods=['GET','POST'])
def search_db():

    """
    Connects user to MySQL server

    Parameters
    ----------
    db: str
        name of ChemBDDB database being queried
    meta: str
        user specifies desired queries using the following keywords:
            'search-query-molecule'
            'search-query-des'
            'next-50'
            'prev-50'
            'all'
            'download_csv'
            'download_json'
            'orderby_property'
        if none of the above keywords are selected, then the function gives blank results
    from_form: dict
        dictionary containing keys and values for each function. Values are lists containing value of below specifications.
        Keys options for each meta selection are as follows:
        'search-query-molecule':
            'property_id' int
            'property_to_val' float
            'property_from_val' float
            'smiles_search' None
            'SMILES' str
            'MW' None
            'MW_to_val' float
            'MW_from_from' float
        'search-query-des':
            'property_id' int
            'property_to_val' float
            'property_from_val' float
            'smiles_search' None
            'HBA_SMILES' str
            'HBA_MW' None
            'HBA_MW_to_val' float
            'HBA_MW_from_val' float
            'HBD_SMILES' str
            'HBD_MW' None
            'HBD_MW_to_val' float
            'HBA_MW_from_val' float
        'next-50':
            No input
        'prev-50':
            No input
        'download_csv' or 'download_json':
            'molecule' None
            'des' None
        'orderby_property':
            'select_order': str
                If value is 'ascending', order will be changed to ascending order        

    Returns
    -------
    sdb: dict
        contains outputs from query and info needed to generate webpage
    """
   
    all_dbs=[]
    sdb=dict()
    cur.execute('show databases;')
    all_dbs_tup=cur.fetchall()
    from_form=request.form
    #print(all_dbs_tup)
    # for i in all_dbs_tup:
    #     if '_chembddb' in i[0] and 'unit_list' not in i[0]:
    #         m=i[0]
    #         all_dbs.append((m[:-9],))
    # db=db[1:-1]
    # db=db+'_chembddb'
    # cur.execute('USE %s'%db)
    # db = db.replace('_chembddb','')
    db=from_form["db"]
    if "_chembddb" not in db:
        db=db+"_chembddb"
    meta=from_form["meta"]
    cur.execute('USE %s'%db)
    cur.execute('Select * from Property')
    properties=cur.fetchall()
    #print(properties)
    cur.execute('Select id,method_name from Model')
    results=cur.fetchall()
    cur.execute("Select * from Functional")
    functionals=cur.fetchall()
    cur.execute("Select * from basis_set")
    basis_sets=cur.fetchall()
    cur.execute("Select * from forcefield")
    forcefields=cur.fetchall()
    methods=[]
    # search_results useful only in case of smiles_search because the same results are used ... single db call
    global search_results
    global sql
    global ini
    global fin
    global n_res
    global noprev
    global nonext
    global keys
    global multiprop
    multiprop = False
    for i in results:
        methods.append(i[1])
    if 'search-query-molecule' in meta:
        keys=[i for i in from_form if '_id' in i]
        min_max_err=False
        min_max_prop=[]
        props=[]
        p = []
        # tuple of tuples to list of tuples
        for pr in properties:
            p.append(list(pr))
        properties = p
        if len(keys)>0:
            sql = 'select value.molecule_id, molecule.SMILES, molecule.MW, Property.Property_str, value.num_value'
            sql = sql + ' from molecule inner join Value on molecule.id=value.molecule_id'
            sql = sql + ' inner join property on property.id=value.property_id'
            sql = sql + ' where '
            counts_q = 'SELECT COUNT(*) '
            counts_q = counts_q + 'from molecule inner join Value on molecule.id=value.molecule_id '
            counts_q = counts_q + ' inner join property on property.id=value.property_id '
            counts_q = counts_q + 'where '
            #print(from_form)
            for k in keys:
                prop_id = int(from_form[k])
                props.append(prop_id)
                from_val = float(from_form[k[:-3]+'_from_val'])
                to_val = float(from_form[k[:-3]+'_to_val'])
                properties[prop_id-1].append(from_val)
                properties[prop_id-1].append(to_val)
                if from_val > to_val:
                    min_max_err=True
                sql = sql[:sql.rfind('where')+6] + 'value.property_id={0} and value.num_value>{1} and value.num_value<{2} and '.format(prop_id,from_val,to_val) + sql[sql.rfind('where')+6:]
                counts_q = counts_q[:counts_q.rfind('where')+6] + 'value.property_id={0} and value.num_value>{1} and value.num_value<{2} and '.format(prop_id,from_val,to_val) + sql[sql.rfind('where')+6:]
                if len(keys)!=0:
                    sql=sql[:-5]
                    counts_q = counts_q[:-5]
                valsid=' and value.property_id in '
                for i in range(len(props)):
                    if i > 0:
                        valsid = valsid + ',' + str(props[i])
                    else:
                        valsid = valsid + '(' + str(props[i])
                    valsid = valsid + ')'
                    sql = sql + valsid
                    counts_q = counts_q+valsid
        else:
            sql = 'SELECT ID, SMILES, MW from Molecule where '
            counts_q = 'SELECT COUNT(*) FROM molecule where '
        MW_to = None
        if 'MW' in from_form:
            from_val=float(from_form['MW_from_val'])
            to_val=float(from_form['MW_to_val'])
            MW_from = from_val
            MW_to = to_val
            if from_val > to_val:
                min_max_err=True
            if len(keys)!=0:
                sql = sql+" and molecule.MW > {} and molecule.MW < {} ".format(float(from_form['MW_from_val']),float(from_form['MW_to_val']))
                counts_q = counts_q + " and molecule.MW > {} and molecule.MW < {} ".format(float(from_form['MW_from_val']),float(from_form['MW_to_val']))
            else:
                sql="select id, SMILES, MW from molecule where molecule.MW > {} and molecule.MW < {} ".format(float(from_form['MW_from_val']),float(from_form['MW_to_val']))
                counts_q = "SELECT COUNT(*) FROM molecule where molecule.MW > {} and molecule.MW < {} ".format(float(from_form['MW_from_val']),float(from_form['MW_to_val']))
                keys.append('MW')
        if 'smiles_search' in from_form:
            if len(keys)==0:
                sql = 'SELECT ID, SMILES FROM molecule '
                sql = sql + 'WHERE SMILES = "%s" '%(from_form['SMILES'])
            
            # insert ff, basis, method, func if statements from chembddb code

            #print(sql)
            #because smiles_search will get all results from db because substructure matching is required
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            counts = -1
        else:
            if 'all' not in meta:
                sql = sql + 'limit 50'
            cur.execute(counts_q)
            counts = cur.fetchall()
            counts = counts[0][0]
        if counts==0:
            sql=sql+';'
            # query creation ends here
            temp_col=[]
            temp_met=[]
            if counts == 0:
                if min_max_err==True:
                    n_res = 'Min value entered is > Max value entered for one of the fields above.'
                    columns = ''
                else:
                    n_res = 'Number of results='+ str(counts)+'. No such candidates exist in your database'
                    columns = ''
                if 'MW' in from_form:
                    sdb = {
                        'MW_from': MW_from,
                        'MW_to': MW_to,
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
                else:
                    sdb = {
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
                print(sdb)
        else:
            # executing the query
            print('\nMySQL search query: ')
            print(sql + '\n')
            cur.execute(sql)
            data1 = cur.fetchall()
            data, columns = post_process(sql, data1,'molecule')
            print(data.head())
            # the following section is just to format the column headers that appear on the html page
            temp_met = []
            temp_col = []
            for c in columns:
                if '-' in c:
                    temp_col.append(c.split('-')[0])
                    if len(c.split('-')[0])>2:
                        temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                    else:
                        temp_met.append(c.split('-')[1])
                else:
                    if 'MW' in c:
                        temp_met.append('pybel')
                    else:
                        temp_met.append('')
                    temp_col.append(c)
            # substructure matching using pybel
            try:
                smi_val = None
                if 'smiles_search' in from_form:
                    smarts = pybel.Smarts(from_form['SMILES'])
                    smi_val = smarts
                    for i in range(len(data)):
                        mol = pybel.readstring("smi",data.loc[i]['SMILES'])
                        smarts.obsmarts.Match(mol.OBMol)
                        if len(smarts.findall(mol))==0:
                            data.drop(i,axis=0,inplace=True)
                    if len(data)==0:
                        n_res='Number of results='+str(len(data))+'.\nNo such candidates exit in your database.'
                    else:
                        n_res=len(data)
                        counts=n_res
                    if n_res>50:
                        search_results=data
                        search_results.columns = data.columns
                        data=data[:51]
                elif len(keys)>1:
                    if len(data)==0:
                        n_res='Number of results= '+str(len(data))+'.\nNo such candidates exist in your database.'
                        print('n_res is ' + n_res)
                    else:
                        n_res = len(data)
                        counts = n_res
                        print('n_res is len(data)')
                    if n_res>50:
                        search_results = data
                        search_results.columns = data.columns
                        data = data[:51]
                        print('n_res > 50')
                    else:
                        search_results = data
                        search_results.columns = data.columns
                else:
                    n_res = counts
                    print('n_res = counts')
            except:
                n_res = 'Invalid SMARTS entered'
                data = pd.DataFrame()
            if len(keys) == 1:
                counts = len(data)
                n_res = counts
            desc=['','']
            columns = []
            # creating tuple of tuples for column headers (required for html page)
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            # calculating statistics for each page
            property_names = []
            for i in properties:
                property_names.append(i[1])
            data = data.convert_dtypes()
            for i in data.columns:
                if i in property_names:
                    desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
            data = tuple(data.itertuples(index=False,name=None))
            if len(columns) == 2:
                to_order = False
            else:
                to_order = True
            noprev=True
            db = db.replace('_chembddb','')
            ini=0
            if type(n_res) != str and n_res <50:
                fin = n_res
                nonext=True
            else:
                nonext=False
                fin=50
            if 'MW' in from_form:
                sdb = {
                    'ini': ini,
                    'fin': fin,
                    'MW_from': MW_from,
                    'MW_to': MW_to,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
            else:
                sdb = {
                    'ini': 0,
                    'fin': fin,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
        print('Molecule search completed!')
        #return(sdb)
    elif 'search-query-des' in meta:
        
        # TODO: allow user to search by smiles or MW of HBA and HBD, ratio, property value
        keys = [i for i in from_form if '_id' in i]
        min_max_err = False
        min_max_prop=[]
        props=[]

        p = []
        # tuple of tuples to list of tuples
        for pr in properties:
            p.append(list(pr))
        properties = p
        if len(keys)>0:
            # TODO: add in basis_set, functionals, model, forcefield
            sql= 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES, MW as HBA_MW FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
            sql = sql + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES, MW as HBD_MW FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
            sql = sql + 'SELECT value.des_id, des.hba_id, des_hba.HBA_SMILES, des_hba.HBA_MW, des.hbd_id, des_hbd.HBD_SMILES,des_hbd.HBD_MW, des.ratio, property.property_str, value.num_value FROM value '
            sql = sql + 'INNER JOIN des on des.id = value.des_id '
            sql = sql + 'INNER JOIN property on property.id = value.property_id '
            sql = sql + 'INNER JOIN des_hba on des_hba.id = value.des_id '
            sql = sql + 'INNER JOIN des_hbd on des_hbd.id = value.des_id '
            sql = sql + 'where '
            counts_q = 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES, MW as HBA_MW FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
            counts_q = counts_q + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES, MW as HBD_MW FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
            counts_q = counts_q + 'SELECT COUNT(*) FROM Value '
            counts_q = counts_q + 'INNER JOIN des on des.id = value.des_id '
            counts_q = counts_q + 'INNER JOIN property on property.id = value.property_id '
            counts_q = counts_q + 'INNER JOIN des_hba on des_hba.id = value.des_id '
            counts_q = counts_q + 'INNER JOIN des_hbd on des_hbd.id = value.des_id '
            counts_q = counts_q + 'where '
            for k in keys:
                prop_id= int(from_form[k])
                props.append(prop_id)
                from_val=float(from_form[k[:-3]+'_from_val'])
                to_val=float(from_form[k[:-3]+'_to_val'])
                #print(prop_id)
                properties[prop_id-1].append(from_val)
                properties[prop_id-1].append(to_val)
                if from_val > to_val:
                    min_max_err=True
                sql = sql[:sql.rfind('where')+6] + 'value.property_id = {0} and value.num_value > {1} and value.num_value < {2} and '.format(prop_id,from_val,to_val)+ sql[sql.rfind('where')+6:]
                counts_q = counts_q[:counts_q.rfind('where')+6] + 'value.property_id = {0} and value.num_value > {1} and value.num_value < {2} and '.format(prop_id,from_val,to_val)+ counts_q[counts_q.rfind('where')+6:]
            if len(keys)!=0:
                sql=sql[:-5]
                counts_q = counts_q[:-5]
            valsid=' and value.property_id in '
            for i in range(len(props)):
                if i > 0:
                    valsid = valsid + ',' + str(props[i])
                else:
                    valsid = valsid + '(' + str(props[i])
            valsid = valsid + ')'
            sql =  sql + valsid + ' '
            counts_q = counts_q + valsid + ';'
        else:
            # TODO: add functionals, model, basis_sets, and forcefield
            sql= 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES, MW as HBA_MW FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
            sql = sql + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES, MW as HBD_MW FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
            sql = sql + 'SELECT des.id, des.hba_id, des_hba.HBA_SMILES, des_hba.HBA_MW, des.hbd_id, des_hbd.HBD_SMILES,des_hbd.HBD_MW, ratio from des '
            sql = sql + 'INNER JOIN des_hba on des_hba.id = des.id '
            sql = sql + 'INNER JOIN des_hbd on des_hbd.id = des.id '
            sql = sql + 'where '
        HBA_MW_to = None
        # sql = 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES, MW as HBA_MW FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
        # sql = sql + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES, MW as HBD_MW FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
        # sql = sql + 'select des.id, des.HBA_id, des_hba.HBA_SMILES, des_hba.HBA_MW, des.HBD_id, des_hbd.HBD_SMILES, des_hbd.HBD_MW, ratio from des inner join des_hba on des_hba.id = des.id '
        # sql = sql + 'inner join des_hbd on des_hbd.id = des.id '
        if 'HBA_MW' in from_form:
            # search for des containing HBA by MW
            HBA_from_val=float(from_form['HBA_MW_from_val'])
            HBA_to_val=float(from_form['HBA_MW_to_val'])
            HBA_MW_from = HBA_from_val
            HBA_MW_to = HBA_to_val
            if HBA_from_val > HBA_to_val:
                min_max_err=True
            if len(keys)!=0:
                sql = sql+" and des_hba.HBA_MW > {} and des_hba.HBA_MW < {} ".format(float(from_form['HBA_MW_from_val']),float(from_form['HBA_MW_to_val']))
            else:
                sql = sql + 'des_hba.HBA_MW > {} and des_hba.HBA_MW < {} '.format(float(from_form['HBA_MW_from_val']),float(from_form['HBA_MW_to_val']))
                keys.append('HBA_MW')
        HBD_MW_to = None
        if 'HBD_MW' in from_form:
            # search for des containing HBD by MW
            HBD_from_val=float(from_form['HBD_MW_from_val'])
            HBD_to_val=float(from_form['HBD_MW_to_val'])
            HBD_MW_from = HBD_from_val
            HBD_MW_to = HBD_to_val
            if HBD_from_val > HBD_to_val:
                min_max_err=True
            if len(keys)!=0:
                sql = sql+" and des_hbd.HBD_MW > {} and des_hbd.HBD_MW < {} ".format(float(from_form['HBD_MW_from_val']),float(from_form['HBD_MW_to_val']))
            else:
                sql = sql + 'des_hbd.HBD_MW > {} and des_hbd.HBD_MW < {} '.format(float(from_form['HBD_MW_from_val']),float(from_form['HBD_MW_to_val']))
                keys.append('HBD_MW')
        if 'smiles_search' in from_form:
            if len(keys)==0:
                sql= 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
                sql = sql + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
                sql = sql + 'SELECT des.id, des.hba_id, des_hba.HBA_SMILES, des.hbd_id, des_hbd.HBD_SMILES, des.ratio from des '
                sql = sql + 'INNER JOIN des_hba on des_hba.id = des.id '
                sql = sql + 'INNER JOIN des_hbd on des_hbd.id = des.id '
                if 'HBA_SMILES' in from_form and 'HBD_SMILES' in from_form:
                    sql = sql + 'WHERE HBA_SMILES = "%s" '%(from_form['HBA_SMILES'])
                    sql = sql + 'AND HBD_SMILES = "%s" '%(from_form['HBD_SMILES'])
                elif 'HBA_SMILES' in from_form and 'HBD_SMILES' not in from_form:
                    sql = sql + 'WHERE HBA_SMILES = "%s" '%(from_form['HBA_SMILES'])
                elif 'HBD_SMILES' in from_form and 'HBA_SMILES' not in from_form:
                    sql = sql + 'WHERE HBD_SMILES = "%s" '%(from_form['HBD_SMILES'])
                
        if 'method' in from_form:
            met_id=0
            for m in results:
                if m[1]==from_form['method_name']:
                    met_id=m[0]
            if len(keys)==0:
                sql=sql+" Value.model_id = {}".format(met_id)
            else:
                sql=sql+" and Value.model_id ={}".format(met_id)
        if 'func' in from_form:
            if len(keys)==0:
                sql=sql+' Value.functional_id={}'.format(from_form['functional_name'])
            else:
                sql=sql+' and Value.functional_id={}'.format(from_form['functional_name'])
        if 'basis' in from_form:
            if len(keys)==0:
                sql=sql+' Value.basis_id={}'.format(from_form['basis_set'])
            else:
                sql=sql+' and Value.basis_id={}'.format(from_form['basis_set'])
        if 'ff' in from_form:
            if len(keys)==0:
                sql=sql+' Value.forcefield_id={}'.format(from_form['forcefield'])
            else:
                sql=sql+' and Value.forcefield_id={}'.format(from_form['forcefield'])
        
        if len(keys) <= 0 or 'HBA_MW' in from_form or 'HBD_MW' in from_form:
            counts_q= 'WITH des_hba AS (SELECT des.id, HBA_id, SMILES as HBA_SMILES, MW as HBA_MW FROM des INNER JOIN molecule on molecule.id = des.hba_id), '
            counts_q = counts_q + 'des_hbd AS (SELECT des.id, HBD_id, SMILES as HBD_SMILES, MW as HBD_MW FROM des INNER JOIN molecule on molecule.id = des.hbd_id) '
            counts_q =counts_q + 'select count(*) '+sql[sql.find('from'):] +';'
        #print('counts_q')
        #print(counts_q)
        # because smiles_search will get all results from the db because substructure matching is required
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            counts = -1
        else:
            if 'all' not in meta:
                sql = sql + 'limit 50'
            # print(sql)
            #print(counts_q)
            cur.execute(counts_q)
            counts = cur.fetchall()
            counts = counts[0][0]

        sql=sql+';'
        print(sql)
        # query creation ends here
        temp_col = []
        temp_met = []
        if counts ==0:
            if min_max_err==True:
                n_res = 'Min value entered is > Max value entered for one of the fields above.'
                columns=''
            else:
                n_res = 'Number of results= '+ str(counts)+'. No such candidates exist in your database.'
                columns=''
            if 'HBA_MW' in from_form and 'HBD_MW' in from_form:
                sdb = {
                        'HBA_MW_from': HBA_MW_from,
                        'HBA_MW_to': HBA_MW_to,
                        'HBD_MW_from': HBD_MW_from,
                        'HBD_MW_to': HBD_MW_to,
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
            elif 'HBA_MW' in from_form and 'HBD_MW' not in from_form:
                sdb = {
                        'HBA_MW_from': HBA_MW_from,
                        'HBA_MW_to': HBA_MW_to,
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
            elif 'HBD_MW' in from_form and 'HBA_MW' not in from_form:
                sdb = {
                        'HBD_MW_from': HBD_MW_from,
                        'HBD_MW_to': HBD_MW_to,
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
            else:
                sdb = {
                        'properties': properties,
                        'columns': columns,
                        'methods': methods,
                        'n_res': n_res,
                        'all_dbs': all_dbs
                    }
        else:
            # executing the query
            print('\nMySQL search query: ')
            print(sql + '\n')
            cur.execute(sql)
            data1 = cur.fetchall()
            data,columns = post_process(sql, data1, 'des')
            print(data.head())
            # the following section is just to format the column headers that appear on the html page
            for c in columns:
                if '-' in c:
                    temp_col.append(c.split('-')[0])
                    if len(c.split('-'))>2:
                        temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                    else:
                        temp_met.append(c.split('-')[1])
                else:
                    if '_MW' in c:
                        temp_met.append('pybel')
                    else:
                        temp_met.append('')
                    temp_col.append(c)
            # substructure matching using pybel
            try:
                smi_val = None
                if 'smiles_search' in from_form:
                    if 'HBA_SMILES' in from_form:
                        # SMARTS only works for the organic compounds, not for the salt. Getting rid of ions to get SMARTS for other part
                        hba = from_form['HBA_SMILES']
                        hba_1, hba_2 = hba.split('.')
                        smartsA = pybel.Smarts(hba_1)
                        #print(smartsA)
                        smi_valA = smartsA
                        for i in range(len(data)):
                            molA = pybel.readstring("smi",data.loc[i]['HBA_SMILES'])
                            smartsA.obsmarts.Match(molA.OBMol)
                            if len(smartsA.findall(molA))==0:
                                data.drop(i,axis=0,inplace=True)
                        if len(data)==0:
                            n_res='Number of results='+ str(len(data))+'\nNo such candidates exist in your database.'
                        else:
                            n_res=len(data)
                            counts = n_res
                        if n_res>50:
                            search_results = data
                            search_results.columns = data.columns
                            data = data[:51]
                        else:
                            search_results = data
                            search_results.columns = data.columns
                    if 'HBD_SMILES' in from_form:
                        smartsD = pybel.Smarts(from_form['HBD_SMILES'])
                        smi_valD = smartsD
                        for i in range(len(data)):
                            molD = pybel.readstring("smi",data.loc[i]['HBD_SMILES'])
                            smartsD.obsmarts.Match(molD.OBMol)
                            if len(smartsD.findall(molD))==0:
                                data.drop(i,axis=0,inplace=True)
                        if len(data)==0:
                            n_res='Number of results='+ str(len(data))+'\nNo such candidates exist in your database.'
                        else:
                            n_res=len(data)
                            counts = n_res
                        if n_res>50:
                            search_results = data
                            search_results.columns = data.columns
                            data = data[:51]
                        else:
                            search_results = data
                            search_results.columns = data.columns
                elif len(keys)>1:
                    if len(data)==0:
                        n_res='Number of results= '+str(len(data))+'.\nNo such candidates exist in your database.'
                        print('n_res is ' + n_res)
                    else:
                        n_res = len(data)
                        counts = n_res
                        print('n_res is len(data)')
                    if n_res>50:
                        search_results = data
                        search_results.columns = data.columns
                        data = data[:51]
                        print('n_res > 50')
                    else:
                        search_results = data
                        search_results.columns = data.columns
                else:
                    n_res = counts
                    print('n_res = counts')
            except:
                n_res = 'Invalid SMARTS entered.'
                #print(n_res)
                data = data
                print('n_res is ' + n_res)
            if len(keys) == 1:
                n_res = counts
            desc = ['','']
            columns = []
            # creating tuple of tuples for column headers (required for html page)
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            # calculating statistics for each page
            property_names = []
            for i in properties:
                property_names.append(i[1])
            data = data.convert_dtypes()
            for i in data.columns:
                if i in property_names:
                    desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
            data = tuple(data.itertuples(index=False,name=None))
            if len(columns) == 2:
                to_order = False
            else:
                to_order = True
            noprev = True
            db = db.replace('_chembddb','')
            ini=0
            if type(n_res) != str and n_res < 50:
                fin = n_res
                nonext=True
            else:
                nonext=False
                fin = 50
            if 'HBA_MW' in from_form and 'HBD_MW' in from_form:
                sdb = {
                    'ini': ini,
                    'fin': fin,
                    'HBA_MW_from': HBA_MW_from,
                    'HBA_MW_to': HBA_MW_to,
                    'HBD_MW_from': HBD_MW_from,
                    'HBD_MW_to': HBD_MW_to,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
            elif 'HBA_MW' in from_form and 'HBD_MW' not in from_form:
                sdb = {
                    'ini': ini,
                    'fin': fin,
                    'HBA_MW_from': HBA_MW_from,
                    'HBA_MW_to': HBA_MW_to,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
            elif 'HBD_MW' in from_form and 'HBA_MW' not in from_form:
                sdb = {
                    'ini': ini,
                    'fin': fin,
                    'HBD_MW_from': HBD_MW_from,
                    'HBD_MW_to': HBD_MW_to,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
            else:
                sdb = {
                    'ini': 0,
                    'fin': fin,
                    'data': data,
                    'properties': properties,
                    'columns': columns,
                    'temp_met': temp_met,
                    'n_res': n_res,
                    'to_order': to_order,
                    'all_dbs': all_dbs,
                    'desc': desc,
                    'noprev': noprev,
                    'nonext': nonext
                }
        print('DES search completed!')
        # return(sdb)
    if 'next-50' in meta:
        nonext=False
        n_res_done = 0
        if ' offset' in sql:
            # checking the offset in the previous sql query, n_res_done tells us how many results have been displayed already
            n_res_done = int(sql[sql.rfind('offset')+7:-1])
            sql = sql[:sql.rfind('offset')+6] + ' '+ str(n_res_done + 50) + ';'
            n_res_done = n_res_done +50
        else:
            n_res_done = 50
            sql = sql[:-1]+' offset ' +str(n_res_done)+';'
        if ('property' not in sql and 'MW' not in sql) or len(keys)>1:
            data = search_results[n_res_done:n_res_done+51]
            columns = list(data.columns)
            if len(keys)==0:
                to_order = False
            else:
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(na/na)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
                columns=[c.replace('(na)','') for c in columns]
                to_order = True
        else:
            to_order = True         
            cur.execute(sql)
            data1=cur.fetchall()
            if 'molecule' in from_form.keys():
                data, columns = post_process(sql, data1,'molecule')
            elif 'des' in from_form.keys():
                data, columns = post_process(sql, data1, 'des')
        temp_col=[]
        temp_met=[]
        for c in columns:
            if '-' in c:
                temp_col.append(c.split('-')[0])
                if len(c.split('-'))>2:
                    temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                else:
                    temp_met.append(c.split('-')[1])
            else:
                temp_col.append(c)
                if 'MW' in c:
                    temp_met.append('pybel')
                else:
                    temp_met.append('')
        
        desc=['','']
        columns =[]
        ini = n_res_done
        if (counts - n_res_done) < 50: 
            fin = counts
            nonext = True
        else:
            nonext = False
            fin = n_res_done + 50                
        noprev=False
        data = data.convert_dtypes()
        if temp_met!=[]:
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            property_names = []
            for i in properties:
                property_names.append(i[1])
            data = data.convert_dtypes()
            for i in data.columns:
                if i in property_names:
                    desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
        if 'order by' in sql:
            if 'DESC' in sql:        
                data=data.sort_values(by=data.columns[-1],ascending=False)
            else:
                data=data.sort_values(by=data.columns[-1])
        data = tuple(data.itertuples(index=False,name=None))
        if 'MW' in sql and 'value.property_id=' not in sql:
            sdb = {
                'ini': ini,
                'fin': fin,
                'to_order': to_order,
                'noprev': False,
                'nonext': nonext,
                'data': data,
                'properties': properties,
                'columns': columns,
                'temp_met': temp_met,
                'methods': methods,
                'n_res': n_res,
                'basis': basis_sets,
                'functionals': functionals,
                'forcefields': forcefields,
                'all_dbs': all_dbs,
                'title': db,
                'desc': desc
            }
        else:
            sdb = {
                'ini': ini,
                'fin': fin,
                'to_order': to_order,
                'noprev': False,
                'nonext': nonext,
                'data': data,
                'properties': properties,
                'columns': columns,
                'temp_met': temp_met,
                'methods': methods,
                'n_res': n_res,
                'basis': basis_sets,
                'functionals': functionals,
                'forcefields': forcefields,
                'all_dbs': all_dbs,
                'title': db,
                'desc': desc
            }
        print('Next 50 available!')
        
    elif 'prev-50' in meta:
        noprev=False
        n_res_done = 0
        if ' offset' in sql:
            # checking the offset in the previous sql query, this number tells us how many results have been displayed already
            n_res_done = int(sql[sql.rfind('offset')+7:-1]) - 50
            sql = sql[:sql.rfind('offset')+6] + ' '+str(n_res_done) +';'
            noprev = False
        if n_res_done ==0:
            noprev = True
            nonext = False
        if('property' not in sql and 'MW' not in sql) or len(keys)>1:
            data = search_results[n_res_done:n_res_done+51]
            columns = list(data.columns)
            if len(keys)==0:
                to_order = False
            else:
                columns=[c.replace('(NA/NA)','') for c in columns]
                columns=[c.replace('(na/na)','') for c in columns]
                columns=[c.replace('(NA)','') for c in columns]
                columns=[c.replace('(na)','') for c in columns]
                to_order = True
        else:
            cur.execute(sql)
            data1=cur.fetchall()
            if 'molecule' in from_form.keys():
                data, columns = post_process(sql, data1,'molecule')
            elif 'des' in from_form.keys():
                data, columns = post_process(sql, data1, 'des')
        temp_col=[]
        temp_met=[]

        for c in columns:
            if '-' in c:
                temp_col.append(c.split('-')[0])
                if len(c.split('-'))>2:
                    temp_met.append(c.split('-')[1]+'-'+c.split('-')[2])
                else:
                    temp_met.append(c.split('-')[1])
            else:
                temp_col.append(c)
                if 'MW' in c:
                    temp_met.append('pybel')
                else:
                    temp_met.append('')
        
        desc = ['','']
        columns = []
        fin = n_res_done + 50
        if temp_met!= []:
            for i in range(len(temp_met)):
                columns.append((temp_col[i],temp_met[i]))
            property_names = []
            for i in properties:
                property_names.append(i[1])
            data = data.convert_dtypes()
            for i in data.columns:
                if i in property_names:
                    desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))
        if 'order by' in sql:
            if 'DESC' in sql:
                data = data.sort_values(by=data.columns[-1],ascending=False)
            else:
                data = data.sort_values(by=data.columns[-1])
        data = tuple(data.itertuples(index=False, name=None))
        ini = n_res_done
        if 'MW' in sql and 'value.property_id=' not in sql:
            sdb = {
                'ini': n_res_done,
                'fin': fin,
                'to_order': to_order,
                'noprev': noprev,
                'nonext': False,
                'data': data,
                'properties': properties,
                'columns': columns,
                'temp_met': temp_met,
                'methods': methods,
                'n_res': n_res,
                'basis': basis_sets,
                'functionals': functionals,
                'forcefields': forcefields,
                'all_dbs': all_dbs,
                'title': db,
                'desc': desc
            }
        else:
            sdb = {
                'ini': n_res_done,
                'fin': fin,
                'to_order': to_order,
                'noprev': noprev,
                'nonext': False,
                'data': data,
                'properties': properties,
                'columns': columns,
                'temp_met': temp_met,
                'methods': methods,
                'n_res': n_res,
                'basis': basis_sets,
                'functionals': functionals,
                'forcefields': forcefields,
                'all_dbs': all_dbs,
                'title': db,
                'desc': desc
            }
        print('Previous 50 available!')
        #return(sdb)

    if 'download_csv' in meta or 'download_json' in meta:
        desc = ['','']
        print('Downloading results...')
        # re-executing query to get all results

        if 'MW' not in sql and 'property' not in sql:
            data = search_results
            columns = search_results.columns
            to_order = False
        else:
            if 'MW' in sql and 'property' not in sql:
                sql = sql[:sql.rindex('limit')]+';'
            else:
                sql = sql[:sql.rindex(')')+1]+';'
            cur.execute(sql)
            all_results = cur.fetchall()
            if 'molecule' in from_form.keys():
                data, columns = post_process(sql, all_results,'molecule')
            elif 'des' in from_form.keys():
                data, columns = post_process(sql,all_results,'des')
            to_order = True

        if 'download_json' in meta:
            import json
            data.to_json('results.json')
            msg = 'Results have been downloaded as results.json!'
        else:
            print(data)
            data.to_csv('results.csv',index=None)
            msg='Results have been downloaded as results.csv!'
        #print(msg)
        # Results that go to the html are still limited to 50
        property_names = []
        for i in properties:
            property_names.append(i[1])
        data = data.convert_dtypes()
        for i in data.columns:
            if i in property_names:
                desc.append('mean={}, std={}, min={}, max={}'.format(data[i].describe()['mean'].round(2),data[i].describe()['std'].round(2),data[i].describe()['min'].round(2),data[i].describe()['max'].round(2)))

        columns =[]
        for i in data.columns:
            if '-' not in i:
                if 'MW' in i:
                    columns.append((i,'pybel'))
                else:
                    columns.append((i,''))
            else:
                if len(i.split('-')) > 2:
                    columns.append((i.split('-')[0],i.split('-')[1]+'-'+i.split('-')[2]))
                else:
                    columns.append((i.split('-')[0],i.split('-')[1]))
        data = tuple(data.itertuples(index=False, name=None))
        sdb = {
            'data': data,
            'ini': ini,
            'fin': fin,
            'to_order': to_order,
            'properties': properties,
            'columns': columns,
            'methods': methods,
            'msg': msg,
            'n_res': n_res,
            'functionals': functionals,
            'basis': basis_sets,
            'forcefields': forcefields,
            'all_dbs': all_dbs,
            'noprev': noprev,
            'nonext': nonext,
            'title': db,
            'desc': desc
        }
        print(msg)
        #return(sdb)
    elif 'orderby_property' in meta:
        ascending = True
        if 'ascending' in from_form['select_order']:
            if 'order by' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('limit')] + ' order by molecule.MW ' +sql[sql.rindex('limit'):]
                else:
                    sql = sql[:sql.rindex(')')+1]+ ' order by Value.num_value ' + sql[sql.rindex(')')+1:]
            elif 'order by' in sql and 'DESC' in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('DESC')] + sql[sql.rindex('DESC')+5:]
                else:
                    sql = sql[:sql.rindex('value')+5] + ' ' + sql[sql.rindex('value')+11:]
        else:
            ascending = False
            if 'order by' in sql and 'DESC' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('MW')+3] + 'DESC ' + sql[sql.rindex('MW')+3:]   
                else:               
                    sql = sql[:sql.rindex('value')+5] + ' DESC' + sql[sql.rindex('value')+5:]
            elif 'order by' not in sql:
                if 'property' not in sql and 'MW' in sql:
                    sql = sql[:sql.rindex('limit')] + ' order by molecule.MW DESC ' +sql[sql.rindex('limit'):]
                else:
                    sql = sql[:sql.rindex(')')+1]+ ' order by value.num_value DESC ' + sql[sql.rindex(')')+1:]
        
        if sql.count('value.num_value') > 4:
            multiprop = True
            sql = sql[:sql.rfind(')')+1]+';'
        cur.execute(sql)
        all_results = cur.fetchall()
        if 'molecule' in from_form.keys():
            data, columns = post_process(sql,all_results,'molecule')
        elif 'des' in from_form.keys():
            data, columns = post_process(sql,all_results,'des')
        search_results = data
        search_results.columns = columns
        if 'MW' not in sql and 'property' in sql:
            search_results = search_results.sort_values(by=from_form['property_orderby'], ascending = ascending)
        desc=['','']
        data = tuple(search_results[:50].itertuples(index=False,name=None))
        for i in search_results.columns[2:]:
            if '-' not in i:
                if 'MW' in i:
                    columns.append((i,'pybel'))
                else:
                    columns.append((i,''))
            else:
                if len(i.split('-')) > 2:
                    columns.append((i.split('-')[0],i.split('-')[1]+'-'+i.split('-')[2]))
                else:
                    columns.append((i.split('-')[0],i.split('-')[1]))
        sdb = {
            'data': data,
            'properties': properties,
            'to_order': True,
            'ini': ini,
            'fin': fin,
            'noprev': noprev,
            'nonext': nonext,
            'columns': columns,
            'methods': methods,
            'n_res': n_res,
            'basis': basis_sets,
            'functionals': functionals,
            'forcefields': forcefields,
            'all_dbs': all_dbs,
            'title': db,
            'desc': desc
        }
        print('Ordered by property!')
        #return(sdb)
    cur.execute('SHOW DATABASES;')
    all_dbs = cur.fetchall()
    return render_template("search.html",sdb=sdb,all_dbs=all_dbs)

@app.route('/delete',methods=['GET','POST'])
def show_databases():
    """
    Shows all databases to user

    Returns
    -------
    all_dbs: list of str
        all databases
    
    """
    cur.execute('SHOW DATABASES;')
    all_dbs = cur.fetchall()
    return(all_dbs)

def drop_database(db):
    """
    Drops any user-specified database

    Parameters
    ----------
    db: str
        name of the database user wishes to drop
    
    Returns
    -------
    all_dbs: list of str
        all remaining databases

    """
    cur.execute('DROP DATABASE `%s`;'%(db))
    cur.execute('SHOW DATABASES;')
    all_dbs = cur.fetchall()
    print('Database %s has been dropped!'%db)
    return(all_dbs)

def drop_table(db,table_name):
    """
    Drops any user-specified table from a user-specified database

    Parameters
    ----------
    db: str
        name of the database user wishes to drop from
    table_name: str
        name of the table user wishes to drop
    
    Returns
    -------
    all_tables: list of str
        list of all tables remaining in database

    """
    cur.execute('DROP TABLE `%s`'%db+'.`%s`;'%(table_name))
    print('Table %s has been dropped!'%(table_name))
    cur.execute('SHOW TABLES FROM `%s`;'%(db))
    all_tables = cur.fetchall()
    return(all_tables)

def drop_row(db,meta,from_form):
    sdb = search_db(db,meta,from_form)
    data = sdb['data']
    val_drop = []
    if 'des' in meta:
        for i in data:
            val_drop.append([i[4],i[1]])
        val_drop = tuple(tuple(x) for x in val_drop)
        for i in val_drop:
            cur.execute('DELETE FROM Value WHERE num_value = %s AND DES_ID = %s;'%i)
    else:
        for i in data:
            val_drop.append([i[3],i[0]])
        val_drop = tuple(tuple(x) for x in val_drop)
        for i in val_drop:
            cur.execute('DELETE FROM Value WHERE num_value = %s AND molecule_id = %s;'%i)
    con.commit()
    print('Selected values have been deleted!')

def backup_restore(host='',user='',pw='',db='',filename=''):
    import subprocess
    if db !='':
        pass
    else:
        try:
            s= 'mysqldump -h localhost -P 3306 -u '+user+' -p'+pw+ ' '+db+' --single-transaction  > ' + filename
            subprocess.Popen(s,shell=True)
            return 'done'
        except:
            return 'error'

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD FOLDER'],filename)

def run_config():
    print('Open http://localhost:5000/connect')
    UPLOAD_FOLDER='./static/files'
    app.config["UPLOAD_FOLDER"]=UPLOAD_FOLDER
    app.run(debug=True,port=5000)
    for i in os.listdir():
        if 'chembddb' in i and 'xyz' in i:
            os.remove(i)
if __name__=="__main__":
    run_config()