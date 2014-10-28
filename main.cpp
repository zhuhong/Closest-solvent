/*
About:
    write all the atoms of the solvent molecule in the range R of the solute to the index file.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <utility>
#include <vector>
#include <map>
#include <iomanip>
#include <stdlib.h>

#include "string_operate.h"
#include "read_ndx.h"
#include "pdb.h"
// #include "my_math.h"

using namespace std;

extern "C"
{
#include "xdrfile_xtc.h"
}

pair<float,int> Min_dist(rvec *x,vector<int> solu_l,float * grid, float * center)
{
    // '''
    // atom_l: A full atom list.
    // solu_l: A solute atom index list.
    // solv_i: Index for one solvent molecule.
    // '''
    float min_dist = 0.0;
    int   min_index =0;

    for (int i=0; i< solu_l.size();i++)
    {

        int item = solu_l[i];
        float tmp = sqrt(pow(x[item-1][0]-grid[0],2) \
                + pow(x[item-1][1]-grid[1],2) \
                + pow(x[item-1][2]-grid[2],2));
        if (i==0)
        {
            min_dist = tmp;
            min_index = item;
            continue;
        }
        if (tmp < min_dist)
        {
            min_dist = tmp;
            min_index = item;
        }  
    }
    pair<float,int> data(min_dist,min_index);

    return data;
}


// def Dist(atom_l,origin,solv_i):
//     dist=(atom_l[solv_i].atom_coor_x-origin[0])**2 \
//     + (atom_l[solv_i].atom_coor_y-origin[1])**2 \
//     + (atom_l[solv_i].atom_coor_z-origin[2])**2

//     return math.sqrt(dist)

// def Get_solvent_list(atom_list):
//     solvent_list = list()
//     for atom in atom_list:
//         if atom_list[atom].residue_name == "WAT" and atom_list[atom].atom_name == "O":
//             solvent_list.append(atom)
//         elif atom_list[atom].residue_name == "SOL" and atom_list[atom].atom_name == "OW":
//             solvent_list.append(atom)

//     return solvent_list

void Get_solute_center(rvec * X_in, vector<int> solute_list, float * center)
{
    for (int i=0;i<3;++)
    {
    	center[i]=0.0;//[0.0,0.0,0.0];
    }
    for (int i=0;i< solute_list.size();i++)
    {
        item = solute_list[i];
        center[0] = center[0] + X_in[item-1][0];
        center[1] = center[1] + X_in[item-1][1];
        center[2] = center[2] + X_in[item-1][2];
    }
    // # print center 
    for (int i=0;i<3;++)
    {
        center[i] = center[i]/solute_list.size();
    }
}


// def Get_Cutoff(atom_list,center,sorted_sol_dist,sol_dict,WAT_NUMBER):
//     CUT_OFF= 0.0
//     for i in range(WAT_NUMBER):
//         temp_id=sol_dict[sorted_sol_dist[i]]
//         dist=((atom_list[temp_id].atom_coor_x-center[0])**2 +\
//             (atom_list[temp_id].atom_coor_y-center[1])**2 +\
//             (atom_list[temp_id].atom_coor_z-center[2])**2)
//         if dist > CUT_OFF:
//             CUT_OFF= dist

//     return math.sqrt(CUT_OFF)*1.2



void pRDF(char * top_file, char * trj_file,char * index_file,char * trjout_file, float CUT_OFF)
{
    /*
    Save the WAT_NUMBER closest water molecules to a new trajectory file.
    */
    // START_TIME       =Time.time()
    
    map<int,atom> atom_list    =read_pdb_to_atom(top_file);
    // atom_list        =copy.deepcopy(Atoms)
    vector <Index_class> index_list       =Read_index_to_Inclass(index_file);
	Print_Index(index_list);
	cout << "Choosing solute:" ;
	cin >> solute_index; 
	vector<int> solute_list  =index_list[solute_index].group_list ;
	cout << "Choosing solvent: " ;
	cin >> solvent_index;
	vector<int> solvent_list = index_list[solvent_index].group_list ;
    int WAT_NUMBER       =0;

    // if len(solvent_list) < WAT_NUMBER:
    //     print "Error: The number of water molecules (%d) is less than the critical number (%d)." \
    //     %(len(solvent_list),WAT_NUMBER)
    //     sys.exit()

	int natoms,step;
	float time_temp;
	float p;
	// vector<double> result;
	matrix box;
	rvec *X_in;
	XDRFILE *xtc_in;
	XDRFILE *xtc_out;
	xtc_in=xdrfile_open(trj_file,"r");
	int read_return=read_xtc_natoms(trj_file,&natoms);
	X_in=(rvec * )calloc(natoms,sizeof(X_in[0]));

	xtc_out = xdrfile_open(trjout_file,"w");



    
    // # OUTPUT_ATOMS =len(solute_atom_list)+3*WAT_NUMBER
    // XTC_write    =libxdrfile.xdrfile_open(trjout_file,'w')
    // # x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32)

// # loop through file until return status signifies end or a problem
// # (it should become exdrENDOFFILE on the last iteration)
    // status      = libxdrfile.exdrOK
    bool FIRST_FRAME =true;
    // TEST_SET    = true;
    // test_count  =0
    // solu_set   =set()
    // OFF     =0.0
    while(1)
    {
        read_return=read_xtc(xtc_in,natoms,&step,&time_temp,box,X_in,&p);
        if(read_return!=0)
		{
			break;
		}
        // # do something with x
        map<float,int> sol_dict;
        vector<float> sol_dist;


        if ( FIRST_FRAME == true)
        {
            for (int i=0;i<solvent_list.size();i++)
            {
            	pair <float,int> min_dist_index; 
            	int solvent = solvent_list[i];
                min_dist_index=Min_dist(X_in,solute_list,solvent);
                float min_dist= min_dist_index.first;
                // min_index = min_dist_index.second;
                sol_dist.push_back(min_dist);
                sol_dict[min_dist]=solvent;
            }
		}
        else
        {
            int _tmp_count =0;
            float center[3];
            Get_solute_center(X_in,solute_list,center);
            for (int i=0;i< solvent_list.size();i++)
            { //solvent in solvent_list:
            	solvent = solvent_list[i];
                float dist = Dist(X_in,center,solvent);
                if (dist  > OFF )
                {
                    continue;
                }

                _dist_2,_solu_2 = Min_dist(X_in,solute_list,solvent);
                sol_dist.append(_dist_2);
                sol_dict[_dist_2]=solvent;
            }
		}

        sorted_sol_dist=sorted(sol_dist);


        if (FIRST_FRAME == True)
        {
            for xx,item in enumerate(sorted_sol_dist):
               if item > CUT_OFF:
                    WAT_NUMBER = xx+1
                    break

            OUTPUT_ATOMS =len(solute_atom_list)+3*WAT_NUMBER
            x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32)
            print "WAT_NUMBER: %d ;" %(WAT_NUMBER)

            SOLUTE_CENTER=Get_solute_center(atom_list,solute_atom_list)
// #           print SOLUTE_CENTER
            OFF = Get_Cutoff(atom_list,SOLUTE_CENTER,sorted_sol_dist,sol_dict,WAT_NUMBER)
		}


        // #########################################
        for i in range(len(solute_atom_list)):
            x_out[i,0] =x[solute_atom_list[i]-1,0]
            x_out[i,1] =x[solute_atom_list[i]-1,1]
            x_out[i,2] =x[solute_atom_list[i]-1,2]    
        for i in range(WAT_NUMBER):
            for j in range(3):
                x_out[len(solute_atom_list)+3*i+j,0] =x[sol_dict[sorted_sol_dist[i]]-1+j,0]
                x_out[len(solute_atom_list)+3*i+j,1] =x[sol_dict[sorted_sol_dist[i]]-1+j,1]
                x_out[len(solute_atom_list)+3*i+j,2] =x[sol_dict[sorted_sol_dist[i]]-1+j,2]

        status2=libxdrfile.write_xtc(XTC_write,step,time,box,x_out,prec)

        if (FIRST_FRAME == true)
        {
            new_list=[atom_list[i] for i in solute_atom_list]
            for i in range(WAT_NUMBER):
                new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]])
                new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]+1])
                new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]+2])

            out_file=string.split(trjout_file,'.')[0]+".pdb"
            Simple_atom.Save_file(out_file,new_list)

            FIRST_FRAME = false;
        }

        NOW_TIME=Time.time()
        BIN_TIME=NOW_TIME-START_TIME
        sys.stderr.write("step %10d, time %10.2f ps, time used: %6.2f s\r" %(step,time,BIN_TIME))
        sys.stderr.flush()
	}
    // # finally close file
    print ""
    libxdrfile.xdrfile_close(XTC)
    // # libxdrfile.xdrfile_close(XTC_write)
}

def For_allsolvent(top_file,trj_file,trjout_file,WAT_NUMBER=500):
    '''
    Save the WAT_NUMBER closest water molecules to a new trajectory file.
    '''
    START_TIME       =Time.time()
    
    Atoms            =Simple_atom.Get_atom_list(top_file)
    atom_list        =copy.deepcopy(Atoms)
    # index_list       =Index.Read_index_to_Inclass(index_file)
    # solute_atom_list =index_list[solute_index].group_list
    solvent_list     =Get_solvent_list(atom_list)

    if len(solvent_list) < WAT_NUMBER:
        print "Error: The number of water molecules (%d) is less than the critical number (%d)." %(len(solvent_list),WAT_NUMBER)

    natoms       = libxdrfile.read_xtc_natoms(trj_file)
    x            = numpy.zeros((natoms, libxdrfile.DIM), dtype=numpy.float32)
    box          = numpy.zeros((libxdrfile.DIM,libxdrfile.DIM), dtype=numpy.float32)
    XTC          = libxdrfile.xdrfile_open(trj_file, 'r')
    
    OUTPUT_ATOMS =3*WAT_NUMBER
    XTC_write    =libxdrfile.xdrfile_open(trjout_file,'w')
    x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32)

# loop through file until return status signifies end or a problem
# (it should become exdrENDOFFILE on the last iteration)
    status      = libxdrfile.exdrOK
    status2     = libxdrfile.exdrOK
    FIRST_FRAME =True
    # center=numpy.array([0.0,0.0,0.0])
    while status == libxdrfile.exdrOK and status2 == libxdrfile.exdrOK:
        status,step,time,prec = libxdrfile.read_xtc(XTC, box, x)
        # do something with x
        sol_dict=dict()
        sol_dist=list()

###############
        for i in atom_list:
            u_atom = MDPackage.Coor.unit_atom.unit_atom(atom_coor_x=x[i-1,0],atom_coor_y=x[i-1,1],atom_coor_z=x[i-1,2])
            atom_list[i]=u_atom

###############
        
        if FIRST_FRAME == True:
            center=numpy.array([box[0,0]/2,box[1,1]/2,box[2,2]/2])
#            FIRST_FRAME = False

        for solvent in solvent_list:
            dist=Dist(atom_list,center,solvent)
            sol_dist.append(dist)
            sol_dict[dist]=solvent


        sorted_sol_dist=sorted(sol_dist)
        STANDARD=sorted_sol_dist[WAT_NUMBER+1]
  
        for i in range(WAT_NUMBER):
            for j in range(3):
                x_out[3*i+j,0] =x[sol_dict[sorted_sol_dist[i]]-1+j,0]
                x_out[3*i+j,1] =x[sol_dict[sorted_sol_dist[i]]-1+j,1]
                x_out[3*i+j,2] =x[sol_dict[sorted_sol_dist[i]]-1+j,2]

        status2=libxdrfile.write_xtc(XTC_write,step,time,box,x_out,prec)

        if FIRST_FRAME == True:
            new_list=list()
            for i in range(WAT_NUMBER):
                new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]])
                new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]+1])
                new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]+2])

            out_file=string.split(trjout_file,'.')[0]+".pdb"
            Simple_atom.Save_file(out_file,new_list)

            FIRST_FRAME = False

        NOW_TIME=Time.time()
        BIN_TIME=NOW_TIME-START_TIME
        if step % 1000 == 0:
            sys.stderr.write("step %10d, time %10.2f ps, time used: %6.2f s\r" %(step,time,BIN_TIME))
            sys.stderr.flush()

    # finally close file
    print ""
    libxdrfile.xdrfile_close(XTC)
    libxdrfile.xdrfile_close(XTC_write)



def Check_args():
    if len(sys.argv) == 6:
        top_file    =sys.argv[1]
        trj_file    =sys.argv[2]
        index_file  =sys.argv[3]
        trjout_file =sys.argv[4]
        # WAT_NUMBER  =int(sys.argv[5])
        CUT_OFF     =float(sys.argv[5])

        index_list=Index.Read_index_to_Inclass(index_file)
        Index.Print_Index(index_list)
        while True:
            try:
                solute_index=int(raw_input("Choose the solute index:"))
                break
            except:
                print "You should input a number."
                continue

        pRDF(top_file,trj_file,index_file,trjout_file,solute_index,CUT_OFF)

    elif len(sys.argv) == 5:
        top_file    =sys.argv[1]
        trj_file    =sys.argv[2]
        trjout_file =sys.argv[3]
        WAT_NUMBER  =int(sys.argv[4])
        For_allsolvent(top_file,trj_file,trjout_file,WAT_NUMBER)

    else: 
        print "Usage: "
        print "For solute,  using: Closest_solvent.py top_file trj_file index_file trjout_file CUT_OFF"
        print "For solvent, using: Closest_solvent.py top_file trj_file trjout_file WAT_NUMBER"



if __name__ == '__main__':
    Check_args()

