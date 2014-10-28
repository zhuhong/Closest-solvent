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
#include <algorithm> // sort
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


pair<float,int> Min_dist(rvec *x,vector<int> solu_l, int solvent)
{
    // '''
    // atom_l: A full atom list.
    // solu_l: A solute atom index list.
    // solv_i: Index for one solvent molecule.
    // '''
    float min_dist = 0.0;
    int   min_index =0;
    float grid[3];
    for(int i=0;i<3;i++)
    {
    	grid[i]= x[solvent-1][i];
    }

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


float Dist(rvec * X_in, float * origin, int solv_i)
{
	float dist ;
	for (int i=0;i<3;i++)
	{
    	dist= pow(X_in[solv_i-1][0] - origin[0],2) \
    	+ pow(X_in[solv_i-1][1] - origin[1],2) \
     	+ pow(X_in[solv_i-1][2] - origin[2],2) ;
	}
    return sqrt(dist);
}
// def Get_solvent_list(atom_list):
//     solvent_list = list()
//     for atom in atom_list:
//         if atom_list[atom].residue_name == "WAT" and atom_list[atom].atom_name == "O":
//             solvent_list.append(atom)
//         elif atom_list[atom].residue_name == "SOL" and atom_list[atom].atom_name == "OW":
//             solvent_list.append(atom)

//     return solvent_list

void Get_solute_center(rvec * X_in, vector<int> solute_list, float center[])
{
    for (int i=0;i<3;i++)
    {
    	center[i]=0.0;//[0.0,0.0,0.0];
    }
    for (int i=0;i< solute_list.size();i++)
    {
        int item = solute_list[i];
        center[0] = center[0] + X_in[item-1][0];
        center[1] = center[1] + X_in[item-1][1];
        center[2] = center[2] + X_in[item-1][2];
    }
    // # print center 
    for (int i=0;i<3;i++)
    {
        center[i] = center[i]/solute_list.size();
    }
}


float Get_Cutoff(rvec * X_in,float * center,vector<float> sol_dist,map<float,int> sol_dict,int WAT_NUMBER)
{
    float CUT_OFF= 0.0;
    for (int i=0; i< WAT_NUMBER; i++)
    {
        int temp_id=sol_dict[sol_dist[i]];

        float dist=(pow(X_in[temp_id-1][0]-center[0] ,2) +\
                    pow(X_in[temp_id-1][1]-center[1], 2) +\
            		pow(X_in[temp_id-1][2]-center[2], 2));
        if (dist > CUT_OFF)
        {
            CUT_OFF= dist;
        }
	}
    return sqrt(CUT_OFF)*1.2;
}


void pRDF(char * top_file, char * trj_file,char * index_file,char * trjout_file, float CUT_OFF)
{
    /*
    Save the WAT_NUMBER closest water molecules to a new trajectory file.
    */
    
    map<int,atom> atom_list    =read_pdb_to_atom(top_file);
    int solute_index;
    int solvent_index;
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
	int OUTPUT_ATOMS;
	float OFF;
	float time_temp;
	float p;
	// vector<double> result;
	matrix box;
	rvec *X_in;
	rvec *X_out;
	XDRFILE *xtc_in;
	XDRFILE *xtc_out;
	xtc_in	= xdrfile_open(trj_file,"r");
	xtc_out = xdrfile_open(trjout_file,"w");
	int read_return=read_xtc_natoms(trj_file,&natoms);
	X_in=(rvec * )calloc(natoms,sizeof(X_in[0]));

	



    
    // # OUTPUT_ATOMS =len(solute_atom_list)+3*WAT_NUMBER
    // XTC_write    =libxdrfile.xdrfile_open(trjout_file,'w')
    // # x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32)

// # loop through file until return status signifies end or a problem
// # (it should become exdrENDOFFILE on the last iteration)
    // status      = libxdrfile.exdrOK
    bool FIRST_FRAME =true;
    // TEST_SET    = true;
    // test_count  =0
    float center[3];
    // OFF     =0.0
    while(1)
    {
        read_return=read_xtc(xtc_in,natoms,&step,&time_temp,box,X_in,&p);
        if(read_return!=0)
		{
			cout << "hello world"<< endl;
			cout << read_return << endl;
			break;
		}
		cout << "reading frame: "<< time_temp<< endl;
		// if (FIRST_FRAME == false)
		// {
		// 	break;
		// }
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
                // cout << min_dist << endl;
                // min_index = min_dist_index.second;
                sol_dist.push_back(min_dist);
                sol_dict[min_dist]=solvent;
            }
		}
        else
        {
            int _tmp_count =0;
            
            Get_solute_center(X_in,solute_list, center);
            for (int i=0;i< solvent_list.size();i++)
            { //solvent in solvent_list:
            	int solvent = solvent_list[i];
                float dist = Dist(X_in,center,solvent);
                if (dist  > OFF )
                {
                    continue;
                }
                pair <float,int> min_dist_index; 
                min_dist_index = Min_dist(X_in,solute_list,solvent);
                float min_dist= min_dist_index.first;
                sol_dist.push_back(min_dist);
                sol_dict[min_dist]=solvent;
            }
		}

		// cout << sol_dist[0]<< endl;
        sort(sol_dist.begin(),sol_dist.end());
        // cout << sol_dist[0]<< endl;
        // cout << sol_dist.size()<< endl;


        if (FIRST_FRAME == true)
        {
            // for xx,item in enumerate(sorted_sol_dist):
            for(int i=0;i< sol_dist.size();++i)
            {
            	float item = sol_dist[i];
                if (item > CUT_OFF)
                {
                    WAT_NUMBER = i+1;
                    break;
                }
            }
            cout << WAT_NUMBER << endl;

            OUTPUT_ATOMS =solute_list.size()+3*WAT_NUMBER;
            cout << OUTPUT_ATOMS << endl;
            // x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32);
            X_out=(rvec * )calloc(OUTPUT_ATOMS,sizeof(X_in[0]));
            // print "WAT_NUMBER: %d ;" %(WAT_NUMBER)

            Get_solute_center(X_in,solute_list,center);
// #           print SOLUTE_CENTER
            OFF = Get_Cutoff(X_in,center,sol_dist,sol_dict,WAT_NUMBER);
		}


        // #########################################
        // for i in range(len(solute_atom_list)):
        for (int i=0;i< solute_list.size();i++)
        {
        	int item = solute_list[i];
            X_out[i][0] =X_in[item-1][0];
            X_out[i][1] =X_in[item-1][1];
            X_out[i][2] =X_in[item-1][2];    
        }
        // for i in range(WAT_NUMBER):
        for(int i=0;i< WAT_NUMBER; i++)
        {
            for (int j=0; j< 3;j++)
            {
            	int solute_size = solute_list.size();
                X_out[solute_size+3*i+j][0] =X_in[sol_dict[sol_dist[i]]-1+j][0];
                X_out[solute_size+3*i+j][1] =X_in[sol_dict[sol_dist[i]]-1+j][1];
                X_out[solute_size+3*i+j][2] =X_in[sol_dict[sol_dist[i]]-1+j][2];
            }
        }

        // status2=libxdrfile.write_xtc(XTC_write,step,time,box,x_out,prec);
        read_return = write_xtc(xtc_out, OUTPUT_ATOMS, step, time_temp, box, X_out, p);
				     // box_xtc, x_xtc, prec_xtc);

                // read_return=read_xtc(xtc_in,natoms,&step,&time_temp,box,X_in,&p);

        if (FIRST_FRAME == true)
        {
            map<int,atom> new_list; 
            for (int i=0;i< solute_list.size(); i++)
			{
				int item = solute_list[i];
				new_list[item] = atom_list[item];
				new_list[item].x = X_in[item-1][0] * 10;
				new_list[item].y = X_in[item-1][1] * 10;
				new_list[item].z = X_in[item-1][2] * 10;
			}
			int solute_size = solute_list.size();
            for (int i=0;i < WAT_NUMBER; i++) 
            {
            	for(int j=0;j<3;j++)
            	{
            		struct atom temp_atom;
            		int serial;
            		serial = sol_dict[sol_dist[i]]+j;
					temp_atom = atom_list[serial];
					temp_atom.x = X_in[serial-1][0] * 10;
					temp_atom.y = X_in[serial-1][1] * 10;
					temp_atom.z = X_in[serial-1][2] * 10;
					new_list[serial] = temp_atom;
                	// new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]])
                	// new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]+1])
                	// new_list.append(atom_list[sol_dict[sorted_sol_dist[i]]+2])
            	}
            }
            string out_file=Split(trjout_file,'.',0)+".pdb";
            write_pdb(new_list,out_file.c_str());

            FIRST_FRAME = false;
        }

        // NOW_TIME=Time.time()
        // BIN_TIME=NOW_TIME-START_TIME
        // sys.stderr.write("step %10d, time %10.2f ps, time used: %6.2f s\r" %(step,time,BIN_TIME))
        // sys.stderr.flush()
	}
    // # finally close file
    // print ""
    xdrfile_close(xtc_in);
    xdrfile_close(xtc_out);
    // # libxdrfile.xdrfile_close(XTC_write)
}

// def For_allsolvent(top_file,trj_file,trjout_file,WAT_NUMBER=500):
//     '''
//     Save the WAT_NUMBER closest water molecules to a new trajectory file.
//     '''
//     START_TIME       =Time.time()
    
//     Atoms            =Simple_atom.Get_atom_list(top_file)
//     atom_list        =copy.deepcopy(Atoms)
//     # index_list       =Index.Read_index_to_Inclass(index_file)
//     # solute_atom_list =index_list[solute_index].group_list
//     solvent_list     =Get_solvent_list(atom_list)

//     if len(solvent_list) < WAT_NUMBER:
//         print "Error: The number of water molecules (%d) is less than the critical number (%d)." %(len(solvent_list),WAT_NUMBER)

//     natoms       = libxdrfile.read_xtc_natoms(trj_file)
//     x            = numpy.zeros((natoms, libxdrfile.DIM), dtype=numpy.float32)
//     box          = numpy.zeros((libxdrfile.DIM,libxdrfile.DIM), dtype=numpy.float32)
//     XTC          = libxdrfile.xdrfile_open(trj_file, 'r')
    
//     OUTPUT_ATOMS =3*WAT_NUMBER
//     XTC_write    =libxdrfile.xdrfile_open(trjout_file,'w')
//     x_out        =numpy.zeros((OUTPUT_ATOMS,libxdrfile.DIM),dtype=numpy.float32)

// # loop through file until return status signifies end or a problem
// # (it should become exdrENDOFFILE on the last iteration)
//     status      = libxdrfile.exdrOK
//     status2     = libxdrfile.exdrOK
//     FIRST_FRAME =True
//     # center=numpy.array([0.0,0.0,0.0])
//     while status == libxdrfile.exdrOK and status2 == libxdrfile.exdrOK:
//         status,step,time,prec = libxdrfile.read_xtc(XTC, box, x)
//         # do something with x
//         sol_dict=dict()
//         sol_dist=list()

// ###############
//         for i in atom_list:
//             u_atom = MDPackage.Coor.unit_atom.unit_atom(atom_coor_x=x[i-1,0],atom_coor_y=x[i-1,1],atom_coor_z=x[i-1,2])
//             atom_list[i]=u_atom

// ###############
        
//         if FIRST_FRAME == True:
//             center=numpy.array([box[0,0]/2,box[1,1]/2,box[2,2]/2])
// #            FIRST_FRAME = False

//         for solvent in solvent_list:
//             dist=Dist(atom_list,center,solvent)
//             sol_dist.append(dist)
//             sol_dict[dist]=solvent


//         sorted_sol_dist=sorted(sol_dist)
//         STANDARD=sorted_sol_dist[WAT_NUMBER+1]
  
//         for i in range(WAT_NUMBER):
//             for j in range(3):
//                 x_out[3*i+j,0] =x[sol_dict[sorted_sol_dist[i]]-1+j,0]
//                 x_out[3*i+j,1] =x[sol_dict[sorted_sol_dist[i]]-1+j,1]
//                 x_out[3*i+j,2] =x[sol_dict[sorted_sol_dist[i]]-1+j,2]

//         status2=libxdrfile.write_xtc(XTC_write,step,time,box,x_out,prec)

//         if FIRST_FRAME == True:
//             new_list=list()
//             for i in range(WAT_NUMBER):
//                 new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]])
//                 new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]+1])
//                 new_list.append(Atoms[sol_dict[sorted_sol_dist[i]]+2])

//             out_file=string.split(trjout_file,'.')[0]+".pdb"
//             Simple_atom.Save_file(out_file,new_list)

//             FIRST_FRAME = False

//         NOW_TIME=Time.time()
//         BIN_TIME=NOW_TIME-START_TIME
//         if step % 1000 == 0:
//             sys.stderr.write("step %10d, time %10.2f ps, time used: %6.2f s\r" %(step,time,BIN_TIME))
//             sys.stderr.flush()

//     # finally close file
//     print ""
//     libxdrfile.xdrfile_close(XTC)
//     libxdrfile.xdrfile_close(XTC_write)




void Print_usage()
{
	cout<<"\t\t"<<":-)  Closest solvent  (-:"<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"This is is a program designed for calculating proximal rdf of solvent."<<endl;
	cout<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout << "Usage: Closest-solvent coor_file traj_file index_file trjout_file CUT_OFF"<< endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"Written by Zhu H. VERSION 1.0.0"<<endl;
	cout<<endl;
}

int main(int argc,char * argv[])
{
	// unit is nm
	char * coor_file;
	char * traj_file;
	char * index_file;
	// char * data_file;
	char * trjout_file;
	float cutoff = 1.0;

	switch(argc)
	{
		case 6:
			coor_file = argv[1];
			traj_file = argv[2];
			index_file = argv[3];
			trjout_file = argv[4];
			cutoff = atof(argv[5]);

			pRDF(coor_file, traj_file, index_file, trjout_file, cutoff);
			// void pRDF(char * top_file, char * trj_file,char * index_file,char * trjout_file, float CUT_OFF)
			break;

		case 2:
			Print_usage();
			exit(0);

		default:
			Print_usage();
			exit(0);
	}
}
