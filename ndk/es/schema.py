from ndk.ui import iprint,wprint,eprint
import ndk.features
# Event Store SQL schemas - for now, this is simply a repository for
# table definitions, and offers callers a convenient way to create
# tables and perform simple selections.  Dictionary-based:

class schema ():
    def __init__(self):
        self.defs = {}
        self.defs['properties']  =  'CREATE TABLE properties ( property varchar(255), value varchar(255) );'
        self.defs['basis']       =  'CREATE TABLE basis ( basisID integer primary key, offset real, sigma real, window real, nfuncs integer, ftype varchar(255) );'
        self.defs['basisname']   =  'CREATE TABLE basisname ( basisID integer primary key, name varchar(255) );'
        self.defs['functions']   =  'CREATE TABLE functions ( basisID integer, component integer, k integer, value real );'
        self.defs['event']       =  'CREATE TABLE event ( eventlabel varchar(255),  samplenum integer);'

        # Spike event tables - note that we have to turn foreign keys
        # ON for the constraints to work:
        self.defs['spiketimes']    =  'CREATE TABLE spiketimes( spikeID integer primary key, samplenum integer );'
        self.defs['spikebasis']    =  'CREATE TABLE spikebasis( spikeID integer primary key, basisID integer, FOREIGN KEY (spikeID) REFERENCES spiketimes(spikeID) ON DELETE CASCADE  );'
        self.defs['spikecoefs']    =  'CREATE TABLE spikecoefs( spikeID integer, channel integer, cindex integer, coef real, FOREIGN KEY (spikeID) REFERENCES spiketimes(spikeID) ON DELETE CASCADE );'
        self.defs['spikechannels'] =  'CREATE TABLE spikechannels( spikeID integer primary key, channel integer );'
        self.defs['spikelabels']   =  'CREATE TABLE spikelabels( spikeID integer primary key, label integer,  FOREIGN KEY (spikeID) REFERENCES spiketimes(spikeID) ON DELETE CASCADE );'

        self.defs['glucose']       =  'CREATE TABLE glucose( samplenum integer primary key, channel integer, value integer );'

        self.defs['meta']          =  'CREATE TABLE meta ( samplerate real, nchannels integer, datafile text);'
        self.defs['history']       =  'CREATE TABLE history( seqno integer primary key, date varchar(255), operation varchar(255) );'
        self.verbose = True

    def create_table(self, es, name):
        es.execute('pragma foreign_key=ON;')
        try:
            es.execute(self.defs[name])
            if self.verbose:
                iprint('Created table {}.'.format(name))
        except Exception as e:
            if self.verbose:
                wprint('Could not create table {}.  Error: {}'.format(name,e))

    # def insert(self, name, values):
    
    # Return a function of three arguments (parameter, window width,
    # nfunctions) that will generate an orthonormal basis set for the

    # specified window.  If basis_generator returns None, then an
    # explicit basis is used and is assumed to be stored in the
    # 'functions' table.
    def basis_generator(self, name):
        if name == 'gaussian_basis':
            return ndk.features.gaussian_basis
        elif name == 'pre_basis':
            return ndk.features.pre_basis
        else:
            return None

