#include "bem.h"

using namespace std;

#define API extern "C" __declspec(dllexport)

API double X0 = -1, X1 = 1, Y0 = -1, Y1 = 1, Z0 = -1, Z1 = 1;
API int nX = 20, nY = 20, nZ = 20;

API D3world* new_world(const char* file, double tol = 0.00001, int maxit = 32, int numMom = 4, int numLev = 6, double spaceunit = 0.001, int segmentation = 1000) {
	return new D3world(file, tol, maxit, numMom, numLev, spaceunit, segmentation);
}

API void del_world(D3world* wr) {
	delete wr;
}

API bool load_dxf(D3world* wr, const char* file, int N) {
	D3ImportedElectrodes *impel = new D3ImportedElectrodes();
	if(!impel->Import(file)) return false;

	for (int i = 0; i < N; ++i) {
		stringstream stm;
		stm << i;
		wr->insert(&(impel->FindElectrode(stm.str().c_str())));
	}

	wr->insert(&(impel->FindElectrode("GROUND")));
	return true;
}

API bool load_csv(D3world* wr, const char* file, char delimeter = ',') {
	ifstream fin(file);
	vector< vector<string> > rows;
	vector<string> row;
	string	line, cell;

	while (getline(fin,line)) {
		stringstream lineStream(line);
		row.clear();

		while( getline( lineStream, cell, ',' ) )
			row.push_back( cell );

		if( !row.empty() )
			rows.push_back( row );
	}
	fin.close();

	/*for( int i=0; i<int(rows.size()); i++ )
	{
		for( int j=0; j<int(rows[i].size()); j++ )
			std::cout << rows[i][j] << " ";

		std::cout << std::endl;
	}*/

	vector< vector<string> >::const_iterator iter = rows.begin();

	int node_size = atoi((*iter)[0].c_str());
	int electrode_size = atoi((*iter)[1].c_str());
	++iter;
	vector<double> X(node_size), Y(node_size), Z(node_size);
	for (int i = 0; i < node_size; ++i, ++iter) {
		X[i] = atof((*iter)[0].c_str());
		Y[i] = atof((*iter)[1].c_str());
		Z[i] = atof((*iter)[2].c_str());
	}
	
	for (int i = 0; i < electrode_size; ++i) {
		int A, B, C;
		D3electrode* el = new D3electrode();
		int element_size = atoi((*iter)[0].c_str());
		if (iter->size() > 1) el->SetVoltage(atof((*iter)[1].c_str()));
		++iter;
		for (int j = 0; j < element_size; ++j, ++iter) {
			A = atoi((*iter)[0].c_str()) - 1;
			B = atoi((*iter)[1].c_str()) - 1;
			C = atoi((*iter)[2].c_str()) - 1;
			el->insert(new D3triangle(X[A],Y[A],Z[A],X[B],Y[B],Z[B],X[C],Y[C],Z[C],true));
		}
		wr->insert(el);
	}
	
	if (iter != rows.end()) {
		X0 = atof((*iter)[0].c_str());
		Y0 = atof((*iter)[1].c_str());
		Z0 = atof((*iter)[2].c_str());
		++iter;
	}
	if (iter != rows.end()) {
		X1 = atof((*iter)[0].c_str());
		Y1 = atof((*iter)[1].c_str());
		Z1 = atof((*iter)[2].c_str());
		++iter;
	} else {
		X1 = X0;
		Y1 = Y0;
		Z1 = Z0;
		X0 = -X0;
		Y0 = -Y0;
		Z0 = -Z0;
	}
	if (iter != rows.end()) {
		nX = atoi((*iter)[0].c_str());
		nY = atoi((*iter)[1].c_str());
		nZ = atoi((*iter)[2].c_str());
	}
	return true;
}

API void solve(D3world* wr) {
	wr->correctNorm(0, 0, 0);
	wr->solve();
}

API void region(D3world* wr, double x0, double x1, int nx, double y0, double y1, int ny, double z0, double z1, int nz) {
	wr->calc(x0,x1,nx,y0,y1,ny,z0,z1,nz);
}

API void region_slow(D3world* wr, double x0, double x1, int nx, double y0, double y1, int ny, double z0, double z1, int nz) {
	wr->calc_slow(x0,x1,nx,y0,y1,ny,z0,z1,nz);
}

API double point(D3world* wr, double x, double y, double z) {
	return wr->calc(x, y, z);
}

API double point_slow(D3world* wr, double x, double y, double z) {
	return wr->calc_slow(x, y, z);
}

#ifdef _WINDLL
API int bem_main(int argc, char* argv[]) {
#else
int main(int argc, char* argv[]) {
#endif
	if (argc < 2) return -1;
	char *path, *base, *ext;
	path = new char[strlen(argv[1]) + 1];
	strcpy(path, argv[1]);
	ext = strrchr(path, '.');
	if (!ext || ext == path) return -1;
	*ext = '\0';
	++ext;
	base = strrchr(path, '\\');
	if (base) {
		if (base != path) ++base;
	} else {
		base = path;
	}
	D3world* wr = new_world(base, 0.000001, 64, 6, 6, 0.001, 1000);
	if (strcmp(ext, "csv") == 0)
		load_csv(wr, argv[1]);
	else if (strcmp(ext, "dxf") == 0) {
		if (argc < 3) return -1;
		load_dxf(wr, argv[1], atoi(argv[2]));
	}
	solve(wr);
	region(wr,X0,X1,nX,Y0,Y1,nY,Z0,Z1,nZ);
	//region_slow(wr,X0,X1,nX,Y0,Y1,nY,Z0,Z1,nZ);
	del_world(wr);
	return 0;
}