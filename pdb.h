#ifndef PDB_H
#define PDB_H

#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>

#include "string_operate.h"


struct  atom
{
	int atom_serial;
	char atom_name[20];
	char residue_name[10];
	int residue_serial;
	float x;
	float y;
	float z;

	char character[10];
	// int atom_serial_number;
	// char atom_name[20];
	char alternate_location_indicator;
	// char reside_name[10];
	char chain_identifier;
	// int reside_sequence_number;
	char code_for_insertions_of_residues;
	// double x_coordinate;
	// double y_coordinate;
	// double z_coordinate;
	double occupancy;
	double temperature_factor;
	char segment_identifier[10];
	char element_symbol[10];
	char charge[10];

};



map<int, atom> read_pdb_to_atom(char * pdb_file_name)
{
	map<int, atom> atom_list;
	int n=0;
	string s;
	ifstream infile(pdb_file_name);
	while(!infile.fail())
	{
		getline(infile,s);
		if(s.find("ATOM")<100||s.find("HETATM")<100)
		{
			struct atom item;
			strcpy(item.character,Split(s.substr(0,6),' ',0).c_str());
			item.atom_serial=atoi(Split(s.substr(6,5),' ',0).c_str());
			strcpy(item.atom_name,Split(s.substr(12,4),' ',0).c_str());
			strcpy(item.residue_name,Split(s.substr(17,3),' ',0).c_str());
			item.residue_serial=atoi(Split(s.substr(22,4),' ',0).c_str());

			atom_list[item.atom_serial]= item;
		}
	}
	infile.close();
	return atom_list;
}

void write_pdb(map<int,atom> atom_list, const char * file_name)
{
	ofstream out(file_name);
	map<int,atom>::iterator it=atom_list.begin();
	for(; it!=atom_list.end(); ++it){
		if(it->second.character!=NULL){
			if(string(it->second.character).size()==4){
				out<<setiosflags(ios::left)<<string(it->second.character)<<"  ";
			}
			else{
				out<<setiosflags(ios::left)<<string(it->second.character);
			}
		}
		else{
			out<<"      ";
		}
		if(it->second.atom_serial!=NULL){
			out<<setiosflags(ios::right)<<setw(5)<<it->second.atom_serial;
		//	cout<<setiosflags(ios::right)<<setw(5)<<it->atom_serial_number;
		}
		else{
			out<<"     ";
		}
		if(it->second.atom_name!=NULL)
		{
			if(string(it->second.atom_name).size()>2){
				out<<setiosflags(ios::left)<<setw(5)<<string(it->second.atom_name);
			}
			if(string(it->second.atom_name).size()==2){
				out<<setiosflags(ios::left)<<setw(4)<<string(it->second.atom_name)<<" ";
			}
			if(string(it->second.atom_name).size()==1){
				out<<setiosflags(ios::left)<<setw(3)<<string(it->second.atom_name)<<"  ";
			}
		//	cout<<setiosflags(ios::left)<<setw(4)<<string(it->atom_name);
		}
		else{
			out<<"    ";
		}
		if(it->second.chain_identifier!=NULL){
			out<<setw(1)<<it->second.alternate_location_indicator;
		//	cout<<setw(1)<<it->alternate_location_indicator;
		}
		else{
			out<<" ";
		}
		if(it->second.residue_name!=NULL){
			out<<setiosflags(ios::right)<<setw(3)<<string(it->second.residue_name);
		//	cout<<setiosflags(ios::right)<<setw(3)<<string(it->reside_name);
		}
		else{
			out<<"   ";
		}
		if(it->second.chain_identifier!=NULL){
			out<<setw(2)<<it->second.chain_identifier;
	//		cout<<setw(2)<<it->chain_identifier;
		}
		else{
			out<<"  ";
		}
		if(it->second.residue_serial!=NULL){
			out<<setiosflags(ios::right)<<setw(4)<<it->second.residue_serial;
		}
		else{
			out<<"    ";
		}
		if(it->second.code_for_insertions_of_residues!=NULL){
			out<<setw(1)<<it->second.code_for_insertions_of_residues;
		}
		else{
			out<<" ";
		}
		if(it->second.x!=NULL){
			out<<fixed<<showpoint;
			out<<setiosflags(ios::right)<<setw(11)<<setprecision(3)<<it->second.x;
		}
		else
		{
			out<<"        ";
		}
		if(it->second.y!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(8)<<setprecision(3)<<it->second.y;
		}
		else
		{
			out<<"        ";
		}
		if(it->second.z!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(8)<<setprecision(3)<<it->second.z;
		}
		else
		{
			out<<"        ";
		}
		if(it->second.occupancy!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(6)<<setprecision(2)<<it->second.occupancy;
		}
		else
		{
			out<<"      ";
		}
		if(it->second.temperature_factor!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(6)<<setprecision(2)<<it->second.temperature_factor;
		}
		else
		{
			out<<"      ";
		}
		if(it->second.segment_identifier!=NULL)
		{
			out<<setiosflags(ios::left)<<setw(11)<<string(it->second.segment_identifier);
		}
		else
		{
			out<<"    ";
		}
		if(it->second.element_symbol!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(2)<<string(it->second.element_symbol);
		}
		else
		{
			out<<"  ";
		}
		if(it->second.charge!=NULL)
		{
			out<<setiosflags(ios::right)<<setw(2)<<string(it->second.charge)<<endl;
		}
		else
		{
			out<<"  "<<endl;;
		}
	}
	out.close();
	cout<<"finish write pdb file "<<file_name<<endl;	
}

#endif


