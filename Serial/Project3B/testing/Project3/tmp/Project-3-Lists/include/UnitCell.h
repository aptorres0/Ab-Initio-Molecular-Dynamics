#ifndef _UNITCELL_H
#define _UNITCELL_H

// Class for n-Dimensional cuboid unit cells
class UnitCell {
public:
  UnitCell() {}
  UnitCell(double _x, double _y, double _z) : lx{_x}, ly{_y}, lz{_z} {}
  UnitCell(Vector3D _l) : lx{_l.GetX()}, ly{_l.GetY()}, lz{_l.GetZ()} {}

private:
  // int nDIM;          // number of dimensions
  double lx, ly, lx; // edge length of each side of the cube
};

#endif