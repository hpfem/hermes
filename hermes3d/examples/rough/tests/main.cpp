#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "../config.h"
#include <hermes3d.h>
#include "../definitions.cpp"

// This test makes sure that the example rough works correctly.

double thresh = 10.e-12;

/////////////////////////////////////////////////////////////////////////////////////////////////////
// dbl: a tampered double for which equality and other binary operators are defined up to a threshold
/////////////////////////////////////////////////////////////////////////////////////////////////////

class dbl
{
	friend std::ostream &operator<<(std::ostream &, const dbl&);
	public:
	double x ;
	dbl(double nx) { x = nx ;} ; 
	dbl() { x=0. ;} ;
	bool operator==(const dbl &rhs) const;
	bool operator!=(const dbl &rhs) const;
	bool operator<(const dbl &rhs) const;
} ;

std::ostream &operator<<(std::ostream &output, const dbl &aaa)
{
	output << aaa.x ;
	return output;
}

bool dbl::operator==(const dbl &rhs) const
{
	if( std::abs(this->x - rhs.x) > thresh ) return false;
	return true;
}

bool dbl::operator!=(const dbl &rhs) const
{
	if( std::abs(this->x - rhs.x) > thresh ) return true;
	return false;
}


bool dbl::operator<(const dbl &rhs) const
{
	if( std::abs(this->x - rhs.x) > thresh )  // they are not equal
	{
		if( this->x < rhs.x ) return(true);
		return(false);
	}
	return false ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// corner is an object that holds the (minimum) x,y coordinates of a column
/////////////////////////////////////////////////////////////////////////////////////////////////////

class corner
{
	friend std::ostream &operator<<(std::ostream &, const corner&);
	public:
	dbl x ;
	dbl y ;

	corner(double nx, double ny) {x = nx ; y = ny ;} ;
	//corner(dbl nx, dbl ny) {x = nx ; y = ny ;} ;
	corner(const corner &copyin);
	~corner(){};
	//corner &operator=(const corner &rhs);
	int operator==(const corner &rhs) const;
	int operator<(const corner &rhs) const;

} ;

corner::corner(const corner &copyin)   // Copy constructor to handle pass by value.
{                            
	x = copyin.x;
	y = copyin.y;
}

std::ostream &operator<<(std::ostream &output, const corner &aaa)
{
	output << "(" << aaa.x << ' ' << aaa.y << ")" << " : " ;
	return output;
}

int corner::operator==(const corner &rhs) const
{
	if( this->x == rhs.x && this->y == rhs.y ) return true;
	return false;
}
 
// This function is required for built-in STL list functions like sort
int corner::operator<(const corner &rhs) const
{
	if( this-> y < rhs.y ) return true ;
	if( this->y == rhs.y && this->x < rhs.x ) return true ;
	return false ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
// we define a map from corners (what identifies columns) to sets of elements associated with them
/////////////////////////////////////////////////////////////////////////////////////////////////////

class column
{
	public:
	double lo ;
	double hi ;
	std::set<unsigned int> elements ;
	column(double nlo, double nhi, std::set<unsigned int> nelements) { this->lo = nlo ; this->hi = nhi ; this->elements = nelements ;} ;
	column(const column &);
	column(); // default constructor
	~column(){};
} ;

column::column(const column &copyin)   // Copy constructor to handle pass by value.
{                            
	this->lo = copyin.lo;
	this->hi = copyin.hi;
	this->elements = copyin.elements ;
}

column::column()
{                            
	this->lo = 0.;
	this->hi = 0.;
}


typedef std::map< corner, column > columnList ; 

// override << opeartor for easy output 
std::ostream &operator<<(std::ostream &output, const std::set<unsigned int> &aaa)
{
	output <<  aaa.size() << " : " ;
	for(std::set<unsigned int>::const_iterator i = aaa.begin() ; i != aaa.end() ; i++) output << *i << " " ;
	return output;
}

int main(int argc, char **args)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh. 
  Mesh mesh;
  H3DReader mloader;
  mloader.load("../lshape_hex.mesh3d", &mesh);

  // rough surface handler
  roughSurface r1 ;
  r1.setFilneName(std::string("1.surf")) ;
  r1.loadSurface() ;
  r1.findMinMax();

  columnList myCols ;
  bool further =true ;
  for(int iter = 0 ; iter < 7 && further == true ; iter++)
  {
    std::cout << "Performing Refinement Level: " << iter << std::endl ;
    further = false ;
    for(std::map<unsigned int, Element*>::iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++)
      if ( it->second->used)
        if (it->second->active)
        {
          std::vector<unsigned int> vtcs(it->second->get_num_vertices()) ;
          //unsigned int vtcs[mesh.elements[idx]->get_num_vertices()] ;
          it->second->get_vertices(&vtcs[0]) ;
          //mesh.vertices[vtcs[0]]->dump() ;
          double minX=+std::numeric_limits<double>::max() ;
          double maxX=-std::numeric_limits<double>::max() ; 
          double minY=+std::numeric_limits<double>::max() ;
          double maxY=-std::numeric_limits<double>::max() ; 
          double minZ=+std::numeric_limits<double>::max() ;
          double maxZ=-std::numeric_limits<double>::max() ; 

          for(unsigned int i=0 ; i<vtcs.size() ; i++)
          {
            double x = mesh.vertices[vtcs[i]]->x ;
            double y = mesh.vertices[vtcs[i]]->y ;
            if(minX > x) minX = x ;
            if(maxX < x) maxX = x ;      
            if(minY > y) minY = y ;
            if(maxY < y) maxY = y ;
          }

          double sX = 0.05e-6;
          double sY = 0.05e-6;
          double deltaX = (maxX-minX) ;
          double deltaY = (maxY-minY) ;
          int nX = (int)floor(deltaX/sX) ;
          int nY = (int)floor(deltaY/sY) ;
          for( int u=0 ; u<=nX ; u++) for( int v=0 ; v<nY ; v++)
          {
            double x = minX + u*sX ;
            double y = minY + v*sY ;
            double z = r1.interpolate(x, y) ;
            if(minZ > z) minZ = z ;
            if(maxZ < z) maxZ = z ;
          }

          myCols[corner(minX,minY)].elements.insert(it->first) ;
          myCols[corner(minX,minY)].lo = minZ ;
          myCols[corner(minX,minY)].hi = maxZ ;

          double deltaZ = std::abs(maxZ - minZ) ;
          if( deltaZ > 2.e-6)
          {
            //std::cout << "\n Iter., element, A.R., deltaZ: " << iter << " " << idx << " " << 2.*3.e-6/deltaX << " " << deltaZ ;
            if(mesh.can_refine_element(it->first, H3D_H3D_REFT_HEX_XY)) mesh.refine_element(it->first, H3D_H3D_REFT_HEX_XY) ;
            further = true ;
          }
          //mesh.elements[idx]->dump() ;
        }  
  }


  for(std::map<unsigned int, Vertex*>::iterator it = mesh.vertices.begin(); it != mesh.vertices.end(); it++)
    if( std::abs(mesh.vertices[it->first]->z - 0.) < 1e-32 ) mesh.vertices[it->first]->z = r1.interpolate(mesh.vertices[it->first]->x,mesh.vertices[it->first]->y) ;	

  for(std::map<unsigned int, Element*>::iterator it = mesh.elements.begin(); it != mesh.elements.end(); it++)
    if ( mesh.elements[it->first]->used) if (mesh.elements[it->first]->active)
    {
      std::vector<unsigned int> vtcs(mesh.elements[it->first]->get_num_vertices()) ;
      mesh.elements[it->first]->get_vertices(&vtcs[0]) ;
      double minZ=+std::numeric_limits<double>::max() ;
      double maxZ=-std::numeric_limits<double>::max() ;
      double deltaZ=-std::numeric_limits<double>::max() ;

      for(std::vector<unsigned int>::iterator i=vtcs.begin() ; i!=vtcs.end() ; i++)
      {
        double z = mesh.vertices[*i]->z ;
        if(minZ > z) minZ = z ;
        if(maxZ < z) maxZ = z ;
        deltaZ = std::abs(maxZ - minZ) ;
      }
      //std::cout << std::endl ;
      if(deltaZ > 4e-6)
        if( mesh.can_refine_element(it->first, H3D_REFT_HEX_Z) ) mesh.refine_element(it->first, H3D_REFT_HEX_Z) ;
    }

  out_mesh_vtk(&mesh, "mesh");

  // missing file 1.surf, supposed to fail temporary.
  return ERR_FAILURE;
}
