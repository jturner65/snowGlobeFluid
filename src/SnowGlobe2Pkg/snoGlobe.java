package SnowGlobe2Pkg;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PShape;

//class to define an instance of a snow globe for rendering and fluid simulation
public class snoGlobe {
	protected static SnowGlobeWin pa;
	protected static mySnowGlobeWin win;
	//object for floor - render to shape 1 time at init, and then just render shape
	public static PShape globeBase;
	//# of cells per dimension
	public static final int snowGlobRad = (int) (SnowGlobeWin.groundRadius*1.05);
	//need to change numDim function - globe radius is much bigger than 10
	public static final int numInDim = 40;
	public static final double cellDim = 1.1 * snowGlobRad / (.5*numInDim);//make bigger than snowglobe rad so bounds are found

	///////////////////
	///simulation parameters
	///////////////////
	
	//vector of e-vex -> array of myVector
	public myVector[] vortVec;//, accelVecX, accelVecY, accelVecZ; used for vort particles
	
	public static final int numCellX = numInDim,numCellY = numInDim, numCellZ = numInDim, 											//corresponds to # of cells in each direction
		numCellXY = numCellX*numCellY, numCells = numCellXY * numCellZ,
		halfNmCellX = (numCellX / 2), halfNmCellY = (numCellY / 2), halfNmCellZ = (numCellZ / 2);
	
	public static final int[] s1i = new int[]{numCellX - 1,numCellY - 1,numCellZ - 1},
							s2i = new int[]{numCellX - 2,numCellY - 2,numCellZ - 2};
	
	//precalced vals to speed up calculations
	public static final double[] 	
		h0sd = new double[] { 1.0 / (2.0* s2i[0]), 1.0 / (2.0* s2i[1]), 1.0 / (2.0* s2i[2])},
		hSd = new double []{ .5 * s2i[0], .5* s2i[1], .5* s2i[2]};
	
	public double diff = .0006,			//set via UI
		visc = 0.0001, 
		deltaT = .01,
		radSq = snowGlobRad * snowGlobRad;			//radSq used for collision if sphere - is square of radius which is 10 in original fluid sim


	public myVector cellSz,						//size of cell in each dimension - defaults to 1
			startLoc,					//location to translate to for start of arrays, to render velocity vectors
			ctrSzHalfNC;				//precalculated (center - sz*halfNumCell)/sz for force query for particles
	public myVectorf		hlfCellSz, center;

	//index in arrays of boundary object, mapped to ptr to fluid bnd obj
	public TreeMap<Integer, myFluidBndObj> sphereBnds;										//cells at boundary of sphere/shaped fluid box

	//use arrays for speed
	public double[] oldDensity, density,Vx, Vy, Vz, Vx0, Vy0, Vz0;
	//whether a cell is oob or not
	public boolean[] isOOB;
	
	///////////////////
	///end simulation parameters
	///////////////////
	
	//constantsdescribing the snow globe
	private static final int numGndTarRings = 12, numFencePosts = 7, postWidth = 7;
		//	numGndSlices = numGndTarRings * numFencePosts;
	private static final float radRat = SnowGlobeWin.groundRadius/(1.0f * numGndTarRings);
	private static final int baseRes = 2048;//resolution of the angles around base, i.e. # of polys.  used to determine how many polys to draw for base
	
	//precalced values specific to the globe base
	private static final float gndRad11 = SnowGlobeWin.groundRadius*1.1f,
			gndRad12 = SnowGlobeWin.groundRadius*1.2f,
			gndRad14 = SnowGlobeWin.groundRadius*1.4f, 
			gndRadOv6 = SnowGlobeWin.groundRadius/6.0f, 
			gndRadOv8 = SnowGlobeWin.groundRadius/8.0f;	

	public snoGlobe(SnowGlobeWin _p, mySnowGlobeWin _win, myVectorf _ctr) {
		pa=_p;
		win = _win;
		initGlobeStand();
		initFluidSim(_ctr);
		initMeshBounds();
	}//ctor
	
	private void initFluidSim(myVectorf _ctr){
		diff = .0006;//control by UI
		visc = 0.0001; 
		deltaT = .01;

		sphereBnds = new TreeMap<Integer, myFluidBndObj>();
		cellSz = new myVector(cellDim,cellDim,cellDim);
		hlfCellSz = new myVectorf(cellDim*.5,cellDim*.5,cellDim*.5);
		initVecs();
		center = new myVectorf(_ctr); 
		setStartLoc();
		myVector tmpHalfCell = new myVector (halfNmCellX, halfNmCellY, halfNmCellZ);
		ctrSzHalfNC= new myVector((center.x - tmpHalfCell.x *cellSz.x)/cellSz.x,
				(center.y - tmpHalfCell.y *cellSz.y)/cellSz.y,
				(center.z - tmpHalfCell.z *cellSz.z)/cellSz.z);			//precalculated (center - sz*halfNumCell)/sz for force query for particles
	}//initFluidSim

	private void initVecs() {
		oldDensity = new double[numCells];
		density = new double[numCells];
		Vx = new double[numCells];
		Vy = new double[numCells];
		Vz = new double[numCells];
		Vx0 = new double[numCells];
		Vy0 = new double[numCells];
		Vz0 = new double[numCells];
		isOOB = new boolean[numCells];
		//init isOOB so that 0 and numCell-1 cells are considered out of bounds
		for(int k = 0; k < numCellZ; ++k){
			for(int j = 0; j < numCellY; ++j){
				for(int i = 0; i < numCellX; ++i){
					if((i==0) || (j==0) || (k==0) || (i==s1i[0]) || (j==s1i[1]) || (k==s1i[2])){isOOB[IX(i,j,k)] = true;}
				}
			}
		}
		
		vortVec = initVecAra(numCells);
//		accelVecX = initVecAra(numCells);
//		accelVecY = initVecAra(numCells);
//		accelVecZ = initVecAra(numCells);
	}
	
	private myVector[] initVecAra(int sz){
		myVector[] vecAra = new myVector[sz];
		for (int i = 0; i < sz; ++i) { vecAra[i] = new myVector(); }
		return vecAra;
	}
	private void setStartLoc() { startLoc = new myVector();startLoc.set(center.x - cellSz.x * (halfNmCellX - .5), center.y - cellSz.y * (halfNmCellY - .5), center.z - cellSz.z * (halfNmCellZ - .5)); }
	///////////////
	//Init Mesh Bounds
	////////////////
//	private myVectorf getDistFromCtrVec(int x, int y, int z) {
//		return new myVectorf(cellSz.x*(x - halfNmCellX), cellSz.y*(y - halfNmCellY), cellSz.z*(z - halfNmCellZ));
//	}

	//x,y,z coords in the world of cell - using halfCell size to find distances in middle of cell
	private myVectorf getIdxLoc(int idx) {
		int[] x = {0}, y = {0}, z = {0};			//x,y,z idx's in 3d version of array
		invIX(idx, x, y, z);
		myVectorf res = new myVectorf(cellSz.x*(x[0] - halfNmCellX), cellSz.y*(y[0] - halfNmCellY), cellSz.z*(z[0] - halfNmCellZ)); //getDistFromCtrVec(x[0], y[0], z[0]);		//res here is coordinate values relative to center (in vector form)
		res._add(center);
		return res;
	}
	//normalized vector from cell @idx to center - length is 1 + difference in distance between min dist vert of cube and threshold/boundary
	private myVectorf getIdxCtrVec(int idx) {return new myVectorf( getIdxLoc(idx),center);}
	
	//if this fluidbox uses a mesh (like a sphere or some other mesh), set up boundary 
	//structure to hold idx of upper fwd left corner of cell containing bound and vector 
	//from center of cell to center of mesh, whose magnitude corresponds 
	//amount of cube inside mesh.
	private void initMeshBounds() {
		double[] sqDistsRad = new double[numCells], sqDistsFloor = new double[numCells];
		double[] amtInside= {0.0};//treating like a pointer so can modify value in checkIDXAtThreshold
		myVectorf distVec = new myVectorf();
		int idx;
		for (int k = 0; k < numCellZ; ++k) {
			for (int j = 0; j < numCellY; ++j) {
				for (int i = 0; i < numCellX; ++i) {
					idx = IX(i, j, k);
					//distVec.set(cellSz.x*(i - halfNmCellX), cellSz.y*(j - halfNmCellY), cellSz.z*(k - halfNmCellZ));			//res here is dist vals from center of cube - assume cube is centered in globe			
					distVec = getIdxCtrVec(idx);
					sqDistsRad[idx] = distVec.sqMagn;// getXYZSqDistFromCtr(i, j, k);
					sqDistsFloor[idx] = distVec.z;		//z component of vector to center
				}
			}
		}
		initIndivMeshBounds(sqDistsRad, radSq );//all distances from the center
		initIndivMeshBounds(sqDistsFloor, center.z+cellSz.z );
	}//initSphereBounds
	
	//initialize the bounds and oob aras with the specific distance arrays
	private void initIndivMeshBounds(double[] sqDists, double thresh) {
		
		double[] amtInside= {0.0};//treating like a pointer so can modify value in checkIDXAtThreshold
		myVectorf distVec = new myVectorf(), loc = new myVectorf();
		int idx;
		for (int k = 1; k < s1i[2]; ++k) {			
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {
					idx = IX(i, j, k);
					amtInside[0] = 0;
					if ((checkIDXAtThreshold(thresh, idx, i, j, k, sqDists, amtInside)) && (!isOOB[idx])) {
//						minA = (minA > amtInside[0] ? amtInside[0] : minA);
//						maxA = (maxA < amtInside[0] ? amtInside[0] : maxA);
						loc.set(getIdxLoc(idx));
						loc._add(hlfCellSz);
						sphereBnds.put(idx,  new myFluidBndObj(idx, getIdxCtrVec(idx)._normalized(), loc , amtInside[0]));			//set idx's in-pointing normal
					} else if(isOOB[idx]) {//remove overlapping bounds from multiple colliders
						sphereBnds.remove(idx);//
					}
				}
			}
		}
		//pa.pr("Fluid Bounds built : min/max amts present inside cubes : " + minA + " | " + maxA);
	}//initSphereBounds
	
	//returns true if cell given by x,y,z encapsulates a given threshold value from the center of the fluid box (i.e. is a boundary cell)
	//sq_dists is an array
	private boolean checkIDXAtThreshold(double sqThresh, int _idx, int _x, int _y, int _z, double[] sq_dists, double[] amtInside) {
		double minVal = 9999999999.0, maxVal = -9999999999.0, d;
		int x, y, z;
		//checking the vertices of each cell - would be better to be half-cell dim
		for (int k = 0; k < 2; ++k) {
			z = _z+k;//(_z >= s1i[2] || _z <= 1 ? _z : _z + k);
			for (int j = 0; j < 2; ++j) {
				y = _y+j;//(_y >= s1i[1] || _y <= 1? _y : _y + j);
				for (int i = 0; i < 2; ++i) {
					x = _x+i;//(_x >= s1i[0] || _x <= 1 ? _x : _x + i);
					d = sq_dists[IX(x, y, z)];
					minVal = (minVal > d ? d : minVal);
					maxVal = (maxVal < d ? d : maxVal);
				}
			}
		}
		amtInside[0] = 0;
		if (minVal >= sqThresh) {//out of bounds for this threshold check
			isOOB[_idx] = true;
			return false;
		} else {//in bounds - check if boundary cell or not - defaults to false, don't set to false since we can then re-call this function for multiple colliders (multiple thresholds and sq_dists)
			if ((maxVal < sqThresh) || (maxVal == minVal)) { 	return false; 	}//not a boundary cell, less than thresh or chk for div by 0 - leave isOOB ara alone
			amtInside[0] = Math.sqrt((sqThresh - minVal)/(maxVal - minVal));//boundary cell
			return true;
		}
	}//checkIDXAtThreshold	
	//for debugging purposes - shows the location in the world and the normal vector for the sphereBnds cells
	public void drawBndsVecs(){
		pa.pushMatrix();pa.pushStyle();
		pa.strokeWeight(.10f / pa.camZScaleMult);
		pa.stroke(255);
		for (Integer key : sphereBnds.keySet()){
			myFluidBndObj bnd = sphereBnds.get(key);
			pa.show(bnd.loc,.25f,SnowGlobeWin.gui_Blue, SnowGlobeWin.gui_Black, false);
//			pa.stroke(255);
//			pa.showVec(bnd.loc, 4.0, bnd.norm);
			pa.stroke(255,0,0,255);
			pa.showVec(bnd.loc, 1.0, bnd.frcAdded);
		}		
		pa.popStyle();pa.popMatrix();
	}

	//will draw the density of the fluid at each cell
	public void drawDensity(){
		pa.pushMatrix();pa.pushStyle();

		double minDens = 8899999999.0, maxDens = -9999999999.0;
		myVectorf loc = new myVectorf();
		pa.noStroke();
		for (int idx = 0; idx < numCells;++idx){
			if (isOOB[idx]) {continue;}
			loc.set(getIdxLoc(idx));
			loc._add(hlfCellSz);			
			pa.drawFluidCell(loc, cellDim, density[idx]);
			minDens = (minDens > density[idx] ? density[idx] : minDens);
			maxDens = (maxDens < density[idx] ? density[idx] : maxDens);
		}
		if((minDens != 0) || (maxDens != 0)){pa.pr("Density ranges : min/max amts present inside cubes : " + minDens + " | " + maxDens);}
		pa.popStyle();pa.popMatrix();
	}//drawDensity
	private static float stM = 200.0f;
	public void drawVelocity(){
		pa.pushMatrix();pa.pushStyle();
		pa.strokeWeight(.10f/ pa.camZScaleMult);
		myVectorf loc = new myVectorf(), loc10 = new myVectorf();// vec = new myVectorf();
		pa.beginShape(PConstants.LINES);
		for (int idx = 0; idx < numCells;++idx){
			if (isOOB[idx]) {continue;}
			loc.set(getIdxLoc(idx));
			loc._add(hlfCellSz);
			loc10.set(loc);
			loc10._add(10.0f*(float)Vx[idx],10.0f*(float)Vy[idx],10.0f*(float)Vz[idx]);
			pa.stroke(0);
			pa.gl_vertex(loc);
			pa.stroke(stM+(float)Vx[idx]*stM,stM+(float)Vy[idx]*stM,stM+(float)Vz[idx]*stM);
			pa.gl_vertex(loc10);
		}
		pa.endShape();
		pa.popStyle();pa.popMatrix();
	}//drawVelocity
	
	public void drawSnowGlobeDebug(){
		if(win.privFlags[mySnowGlobeWin.showBndsIDX]){drawBndsVecs();}
		if(win.privFlags[mySnowGlobeWin.showDensIDX]){drawDensity();}
		if(win.privFlags[mySnowGlobeWin.showVelIDX]){drawVelocity();}
	}
	
	public void drawSnowGlobe(){
		pa.shape(globeBase);    //draw snowglobe base object showBndsIDX,showDensIDX,showVelIDX
		drawGlobe();			//needs to be last or will blot out snowmen	
	}
	
	private double[] tmpDAra;
	//execute a single time step of the fluid simulation
	public void timeStep(){
		//add motion sources into arrays (motion added into V*0/oldDensity)
		//addAllSources();
		tmpDAra = Vx; Vx=Vx0;Vx0=tmpDAra;tmpDAra = Vy; Vy=Vy0;Vy0=tmpDAra;tmpDAra = Vz; Vz=Vz0;Vz0=tmpDAra;
		diffuseSphere(Vx, Vx0, Vy, Vy0, Vz, Vz0, visc, 8);
		projectSphere(Vx, Vy, Vz, Vx0, Vy0, 8);
		
		tmpDAra = Vx; Vx=Vx0;Vx0=tmpDAra;tmpDAra = Vy; Vy=Vy0;Vy0=tmpDAra;tmpDAra = Vz; Vz=Vz0;Vz0=tmpDAra;
		advectSphere(Vx, Vx0, Vy, Vy0, Vz, Vz0);
		projectSphere(Vx, Vy, Vz, Vx0, Vy0, 8);
		
		tmpDAra = density; density=oldDensity;oldDensity=tmpDAra;
		diffSphDens(density, oldDensity, diff, 8, numCellX);
		tmpDAra = density; density=oldDensity;oldDensity=tmpDAra;
		advSphDens(density, oldDensity, Vx, Vy, Vz);
		//vort confine here
		vorticityConfinement(oldDensity);
		
		resetOldVals();//clear out old values so that any UI forces can be added to Vx0, etc
	}//timeStep


	private int forceIDXBnd(int idx, int idxMax, int idxMin) { return (((idx > idxMin)&(idx < idxMax) ? idx : (idx > idxMin) ? idxMax : idxMin)); }
	private double forceIDXBndD(double val, double ubnd, double lbnd) { return (((val > lbnd) & (val < ubnd) ? val : (val  > lbnd) ? ubnd : lbnd)); }
	//private void addAllSources(){for (int i = 0; i < numCells; ++i) { if (isOOB[i]) { continue; } Vx[i] += Vx0[i];Vy[i] += Vy0[i];Vz[i] += Vz0[i];density[i] += oldDensity[i]; }	}
	private int IX(int x, int y, int z) {return x + (y * numCellX) + (z * numCellXY);	}
	//given an index, return the x,y,z coords
	private void invIX(int idx, int[] x, int[] y, int[] z) {
		z[0] = (int)(idx / numCellXY);
		int tmp = (idx % numCellXY);
		y[0] = (int)(tmp/numCellX); 
		x[0] = (int)(tmp % numCellX);
	}

	//solve for all dirs simulataneously
	private void lin_solveSphere(double[] x, double[] x0, double[] y, double[] y0, double[] z, double[] z0, double a, double c, int iter) {
		if (a == 0) {
			System.arraycopy(x0, 0, x, 0, x0.length);
			System.arraycopy(y0, 0, y, 0, y0.length);
			System.arraycopy(z0, 0, z, 0, z0.length);
			set_bndSphere3(x, y, z);
		} else {
			int idx0, idx1, idx2, idx3, idx4, idx5, idx6;
			for (int itr = 0; itr < iter; ++itr) {
				for (int k = 1; k < s1i[2]; ++k) {
					for (int j = 1; j < s1i[1]; ++j) {
						for (int i = 1; i < s1i[0]; ++i) {
							idx0 = IX(i, j, k);
							if (isOOB[idx0]) { continue; }
							idx1 = IX(i + 1, j, k);
							idx2 = IX(i - 1, j, k);
							idx3 = IX(i, j + 1, k);
							idx4 = IX(i, j - 1, k);
							idx5 = IX(i, j, k + 1);
							idx6 = IX(i, j, k - 1);
							x[idx0] = (x0[idx0] + a * (x[idx1] + x[idx2] + x[idx3] + x[idx4] + x[idx5] + x[idx6])) / c;
							y[idx0] = (y0[idx0] + a * (y[idx1] + y[idx2] + y[idx3] + y[idx4] + y[idx5] + y[idx6])) / c;
							z[idx0] = (z0[idx0] + a * (z[idx1] + z[idx2] + z[idx3] + z[idx4] + z[idx5] + z[idx6])) / c;
						}//for i
					}//for j
				}//for k
				set_bndSphere3(x, y, z);
			}//for itr
		}//if a != 0
	}//lin_solveSphere


	private void diffuseSphere(double[] x, double[] x0, double[] y, double[] y0, double[] z, double[] z0, double viscdiff, int iter){
		//void myFluidBox::diffuseSphere(double[] x, double[] x0, double[] y, double[] y0, double[] z, double[] z0, double viscdiff, int iter) {
		double delVisc = (deltaT*viscdiff*(numCells));
		lin_solveSphere(x, x0, y, y0, z, z0, delVisc, 1 + 6 * delVisc, iter);
	}
	
	//handle advection for passed arrays of velocities
	private void advectSphere(double[] _velx, double[] _velx0, double[] _vely, double[] _vely0, double[] _velz, double[] _velz0) {
		double s0, s1, t0, t1, u0, u1,x, y, z;

		int i0, i1, j0, j1, k0, k1;
		//precalced idx's
		int idx0, idx1, idx2, idx3, idx4, idx5, idx6, idx7;

		double dtx = deltaT * (s2i[0]),
			dty = deltaT * (s2i[1]),
			dtz = deltaT * (s2i[2]);
		int IXidx;

		for (int k = 1; k < s1i[2]; ++k) {
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {
					IXidx = IX(i, j, k);
					if (isOOB[IXidx]) { continue; }

					x = i - dtx * _velx0[IXidx];
					x = forceIDXBndD(x, s2i[0] + 0.5, .5);
					i0 = (int)Math.floor(x);	i1 = i0 + 1;  s1 = x - i0;	s0 = 1.0 - s1;

					y = j - dty * _vely0[IXidx];
					y = forceIDXBndD(y, s2i[1] + 0.5, .5);
					j0 = (int)Math.floor(y);	j1 = j0 + 1;	t1 = y - j0;	t0 = 1.0 - t1;

					z = k - dtz * _velz0[IXidx];
					z = forceIDXBndD(z, s2i[2] + 0.5, .5);
					k0 = (int)Math.floor(z);	k1 = k0 + 1;	u1 = z - k0;	u0 = 1.0 - u1;

					idx0 = IX(i0, j0, k0);idx1 = IX(i0, j0, k1);idx2 = IX(i0, j1, k0);idx3 = IX(i0, j1, k1);
					idx4 = IX(i1, j0, k0);idx5 = IX(i1, j0, k1);idx6 = IX(i1, j1, k0);idx7 = IX(i1, j1, k1);

					_velx[IXidx] =	
						s0 * (t0 * (u0 * _velx0[idx0] + u1 * _velx0[idx1]) + (t1 * (u0 * _velx0[idx2] + u1 * _velx0[idx3]))) +
						s1 * (t0 * (u0 * _velx0[idx4] + u1 * _velx0[idx5]) + (t1 * (u0 * _velx0[idx6] + u1 * _velx0[idx7])));

					_vely[IXidx] =
						s0 * (t0 * (u0 * _vely0[idx0] + u1 * _vely0[idx1]) + (t1 * (u0 * _vely0[idx2] + u1 * _vely0[idx3]))) +
						s1 * (t0 * (u0 * _vely0[idx4] + u1 * _vely0[idx5]) + (t1 * (u0 * _vely0[idx6] + u1 * _vely0[idx7])));

					_velz[IXidx] =
						s0 * (t0 * (u0 * _velz0[idx0] + u1 * _velz0[idx1]) + (t1 * (u0 * _velz0[idx2] + u1 * _velz0[idx3]))) +
						s1 * (t0 * (u0 * _velz0[idx4] + u1 * _velz0[idx5]) + (t1 * (u0 * _velz0[idx6] + u1 * _velz0[idx7])));
				}//for i
			}//for j
		}//for k
		set_bndSphere3(_velx, _vely, _velz);
	}//advectSphere

	//diffuse density - 1 d but use sphere bnds
	private void diffSphDens(double[] x, double[] x0, double viscdiff, int iter, int _numCells) {
		double a = (deltaT*viscdiff*(numCells)), c = (1 + 6 * a);
		int idx;
		if(a==0){
			System.arraycopy(x0, 0, x, 0, x0.length);
			set_bndDiffSphere(x);
		} else {
			for (int itr = 0; itr < iter; ++itr) {
				for (int k = 1; k < s1i[2]; ++k) {
					for (int j = 1; j < s1i[1]; ++j) {
						for (int i = 1; i < s1i[0]; ++i) {
							idx = IX(i, j, k);
							if (isOOB[idx]) { continue; }
							x[idx] = (x0[idx] + a *
								(x[IX(i + 1, j, k)]
									+ x[IX(i - 1, j, k)]
									+ x[IX(i, j + 1, k)]
									+ x[IX(i, j - 1, k)]
									+ x[IX(i, j, k + 1)]
									+ x[IX(i, j, k - 1)])) / c;
						}//for i
					}//for j
				}//for k
				set_bndDiffSphere(x);
			}
		}//for itr
	}//diffSphDens
	
	
	//advect density through sphere, using sphere bounds
	private void advSphDens(double[] d, double[] d0, double[] velocX, double[] velocY, double[] velocZ) {
		double s0, s1, t0, t1, u0, u1, x, y, z, dtx = deltaT * (s2i[0]),dty = deltaT * (s2i[1]),dtz = deltaT * (s2i[2]);
		int i0, i1, j0, j1, k0, k1, i, j, k, IXidx;

		for (k = 1; k < s1i[2]; ++k) {
			for (j = 1; j < s1i[1]; ++j) {
				for (i = 1; i < s1i[0]; ++i) {
					IXidx = IX(i, j, k);
					if (isOOB[IXidx]) { continue; }
					//advection grid
					x = i - dtx * velocX[IXidx];
					y = j - dty * velocY[IXidx];
					z = k - dtz * velocZ[IXidx];
					x = forceIDXBndD(x, s2i[0] + 0.5, .5);
					i0 = (int)Math.floor(x);	i1 = i0 + 1;  s1 = x - i0;	s0 = 1.0 - s1;

					y = forceIDXBndD(y, s2i[1] + 0.5, .5);
					j0 = (int)Math.floor(y);	j1 = j0 + 1;	t1 = y - j0;	t0 = 1.0 - t1;

					z = forceIDXBndD(z, s2i[2] + 0.5, .5);
					k0 = (int)Math.floor(z);	k1 = k0 + 1;	u1 = z - k0;	u0 = 1.0 - u1;

					d[IXidx] =
						s0 * (t0 * (u0 * d0[IX(i0, j0, k0)]	+ u1 * d0[IX(i0, j0, k1)]) + (t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)]))) + 
						s1 * (t0 * (u0 * d0[IX(i1, j0, k0)]	+ u1 * d0[IX(i1, j0, k1)]) + (t1 * (u0 * d0[IX(i1, j1, k0)]	+ u1 * d0[IX(i1, j1, k1)])));
				}//for i
			}//for j
		}//for k
		set_bndDiffSphere(d);
	}//advSphDens
	
	private void projectSphere(double[] velocX, double[] velocY, double[] velocZ, double[] p, double[] div, int iter) {
		int idx;
		for (int k = 1; k < s1i[2]; ++k) {
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {
					idx = IX(i, j, k);
					p[idx] = 0;
					if (isOOB[idx]) {	div[idx] = 0;	continue; }
					div[idx] = -(
						(velocX[IX(i + 1, j, k)] - velocX[IX(i - 1, j, k)]) * h0sd[0]
						+(velocY[IX(i, j + 1, k)] - velocY[IX(i, j - 1, k)]) * h0sd[1]
						+(velocZ[IX(i, j, k + 1)] - velocZ[IX(i, j, k - 1)]) * h0sd[2]
						);
				}
			}
		}
		set_bndDiffSphere(div);
		set_bndDiffSphere(p);
		//integration - gauss seidel
		for (int itr = 0; itr < iter; ++itr) {
			for (int k = 1; k < s1i[2]; ++k) {
				for (int j = 1; j < s1i[1]; ++j) {
					for (int i = 1; i < s1i[0]; ++i) {
						idx = IX(i, j, k);
						if (isOOB[idx]) { continue; }
						p[idx] = (div[idx] + (p[IX(i + 1, j, k)] + p[IX(i - 1, j, k)] + p[IX(i, j + 1, k)] + p[IX(i, j - 1, k)] + p[IX(i, j, k + 1)] + p[IX(i, j, k - 1)])) / 6.0;
					}//for i
				}//for j
			}//for k
			set_bndDiffSphere(p);
		}//for each iteration
		
		for (int k = 1; k < s1i[2]; ++k) {
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {
					idx = IX(i, j, k);
					if (isOOB[idx]) { continue; }
					velocX[idx] -= hSd[0] *(p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
					velocY[idx] -= hSd[1] *(p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
					velocZ[idx] -= hSd[2] *(p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
				}
			}
		}
		set_bndSphere3(velocX, velocY, velocZ);
	}//projectSphere

	//vorticity confinement - add back vorticity details lost through numerical dissipation
	//vortN is unused array to hold calcs
	private void vorticityConfinement(double[] vortN) {
		int idx, idx_ijp1k, idx_ijm1k, idx_ip1jk, idx_im1jk, idx_ijkm1, idx_ijkp1;
		double vortEps = deltaT * .01;	//TODO change to allow for user input
		for (int k = 1; k < s1i[2]; ++k) {
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {
					idx = IX(i, j, k);
					if (isOOB[idx]) { continue; }
					idx_ip1jk = IX(i + 1, j, k); idx_im1jk = IX(i - 1, j, k); idx_ijp1k = IX(i, j + 1, k); idx_ijm1k = IX(i, j - 1, k); idx_ijkp1 = IX(i, j, k + 1); idx_ijkm1 = IX(i, j, k - 1);
					//curl operation del cross u -> partial z w/respect to y is the finite diff of the z vels across the y coords
					vortVec[idx].set( 
						((Vy[idx_ijkp1] - Vy[idx_ijkm1]) * h0sd[1]) - ((Vz[idx_ijp1k] - Vz[idx_ijm1k]) * h0sd[2]),
						((Vz[idx_ip1jk] - Vz[idx_im1jk]) * h0sd[2]) - ((Vx[idx_ijkp1] - Vx[idx_ijkm1]) * h0sd[0]),
						((Vx[idx_ijp1k] - Vx[idx_ijm1k]) * h0sd[0]) - ((Vy[idx_ip1jk] - Vy[idx_im1jk]) * h0sd[1]));
				}
			}
		}
		set_bndDiffVec(vortVec);
		for (idx = 0; idx < numCells; ++idx) { vortN[idx] = (isOOB[idx]) ?  0 : vortVec[idx].magn; }

		myVector eta = new myVector(), vf = new myVector();
		
		for (int k = 1; k < s1i[2]; ++k) {
			for (int j = 1; j < s1i[1]; ++j) {
				for (int i = 1; i < s1i[0]; ++i) {	
					idx = IX(i, j, k);
					if (vortN[idx] < .0000001) {	continue;}
					vortVec[idx]._normalize();
					eta.set( ((vortN[IX(i + 1, j, k)] - vortN[IX(i - 1, j, k)]) * h0sd[0]), ((vortN[IX(i, j + 1, k)] - vortN[IX(i, j - 1, k)]) * h0sd[1]), ( (vortN[IX(i, j, k + 1)] - vortN[IX(i, j, k - 1)]) * h0sd[2]));
					eta._normalize();
					vf.set(eta._cross(vortVec[idx]));
					vf._mult(vortEps);
					//cout << "Vx " << idx << " before :  " << Vx[idx] << " vortN : " << vortN[idx] << " invDivX : " << invDivX << " eta : " << eta(0) << "," << eta(1) << "," << eta(2) << " vf : " << vf(0) << "," << vf(1) << "," << vf(2) << " vort : " << vort(0) << "," << vort(1) << "," << vort(2);
					Vx[idx] += vf.x * s2i[0];
					Vy[idx] += vf.y * s2i[1];
					Vz[idx] += vf.z * s2i[2];	
					//cout << "Vx " << idx << " after :  " << Vx[idx]<<endl;
				}
			}
		}
		set_bndSphere3(Vx, Vy, Vz);
	}//vorticityConfinement
	
	private void set_bndDiffSphere(double[] x) {for (Integer key : sphereBnds.keySet()) { x[key] *= sphereBnds.get(key).mag; }}//scale to amt of cube in bounds	
	//private void set_bndDiffSphere(double[] x) {for (Integer key : sphereBnds.keySet()) { x[key] *= sphereBnds.get(key).oneMMag; }}//scale to amt of cube in bounds - boundary conditions, decreases more across larger cubes
	private void set_bndDiffVec(myVector[] egVec) {
		myVector velNorm = new myVector(0, 0, 0);	
		double dotProd;
		myFluidBndObj bnd ;
		for (Integer key : sphereBnds.keySet()) {
			bnd = sphereBnds.get(key);
			dotProd = egVec[key]._dot(bnd.norm);
			if(dotProd >= 0){
//				egVec[key].set(0,0,0);
				continue;
			}//vel in dir of normal
//			velNorm = myVector._mult(bnd.norm,(dotProd * bnd.mag));			//velocity in the normal direction toward center of sphere
			dotProd *= bnd.mag; ;
			velNorm = myVector._mult(bnd.norm,dotProd );			//velocity in the normal direction toward center of sphere
			egVec[key]._sub(velNorm);
			egVec[key]._mult(dotProd);
		}
	}
	
	//address boundary layer values
	private void set_bndSphere3(double[] x, double[] y, double[] z) {
		myVector velVec = new myVector(0, 0, 0), velNorm = new myVector(0, 0, 0);
		myFluidBndObj bnd;
		double dotProd;
		for (Integer key : sphereBnds.keySet()) {
			bnd = sphereBnds.get(key);
			velVec.set( x[key], y[key], z[key]);
			dotProd = velVec._dot(bnd.norm);
			if(dotProd >= 0){
				
				continue;}//vel in dir of normal
//			velNorm = myVector._mult(bnd.norm,(dotProd * bnd.mag));// (velVec.dot(it->second->norm)* it->second->mag)  *it->second->norm;			//velocity in the normal direction toward center of sphere
			dotProd *= bnd.mag; ;
			velNorm = myVector._mult(bnd.norm,dotProd);// (velVec.dot(it->second->norm)* it->second->mag)  *it->second->norm;			//velocity in the normal direction toward center of sphere
			x[key] -=  velNorm.x;
			y[key] -=  velNorm.y;
			z[key] -=  velNorm.z;
			x[key] *=  dotProd;
			y[key] *=  dotProd;
			z[key] *=  dotProd;
		}
	}//set_bndSphere

	public void myFluidBoxAddDensity(int x, int y, int z, double amount) { density[IX(x, y, z)] += amount; }

	public void resetOldVals() {
		pa.dAraFill(Vx0, 0);
		pa.dAraFill(Vy0, 0);
		pa.dAraFill(Vz0, 0);
		pa.dAraFill(oldDensity, 0);
		for (Integer key : sphereBnds.keySet())	{sphereBnds.get(key).frcAdded.set(0,0,0);}
	}
	
	//add the force described by the vector amt to all the boundary cells in the globe, dotted against the normal
	public void myFluidBoxAddDens(double amt, myPoint cellLoc){
		//add density at specified location
		int idx = IX(forceIDXBnd((int)(cellLoc.x), s1i[0],0),
				 forceIDXBnd((int)(cellLoc.y), s1i[1],0),
				 forceIDXBnd((int)(cellLoc.z), s1i[2],0));
		if(isOOB[idx]){return;}
		density[idx] += amt;
	
	}//myFluidBoxAddDens
	
	//add the force described by the vector amt to all the boundary cells in the globe, dotted against the normal
	public void myFluidBoxAddForce(myVector amt){
		myFluidBndObj bnd;
		//myVector amtNorm = myVector._normalize(amt);
		//sphereBnds
		for (Integer key : sphereBnds.keySet()) {
			bnd = sphereBnds.get(key);
			double dotProd = amt._dot(bnd.norm);
			if(dotProd <= 0){bnd.frcAdded.set(0,0,0);continue;}
			bnd.frcAdded.set(myVectorf._mult(bnd.norm, dotProd));
//			//pa.pr("before IDX : " + key + " norm : " + bnd.norm.toStrBrf()+"\tamt:"+amt.toString() + "\tdot prod : " + dotProd + "\t x:"+Vx[key] + " y:"+Vy[key] + " z:"+Vz[key] );
			Vx[key] += bnd.frcAdded.x;
			Vy[key] += bnd.frcAdded.y;
			Vz[key] += bnd.frcAdded.z;
			
//			//try subtraction since setBnds flips
//			Vx[key] -= bnd.frcAdded.x;
//			Vy[key] -= bnd.frcAdded.y;
//			Vz[key] -= bnd.frcAdded.z;
			
//			Vx0[key] = amt.x * dotProd;
//			Vy0[key] = amt.y * dotProd;
//			Vz0[key] = amt.z * dotProd;
			//pa.pr("after x:"+Vx[key] + " y:"+Vy[key] + " z:"+Vz[key] );
		}		
	}//myFluidBoxAddForce

	public void myFluidBoxAddForce(myVector cellLoc, myVector amount) {
		//pa.pr("force addition location in cube : "+cellLoc.toStrBrf());
		int idx = IX(forceIDXBnd((int)(cellLoc.x), s1i[0],0),
					 forceIDXBnd((int)(cellLoc.y), s1i[1],0),
					 forceIDXBnd((int)(cellLoc.z), s1i[2],0));
		if(isOOB[idx]){return;}
		Vx0[idx] = amount.x;
		Vy0[idx] = amount.y;
		Vz0[idx] = amount.z;
	}

	//converts world location into cell locations (int val being cell idx in that axis, float val being interpolant of cell dim along that axis
	public myVector transLocIntoCell(myVector testLoc){return new myVector(ctrSzHalfNC,new myVector (testLoc.x/cellSz.x, testLoc.y/cellSz.y, testLoc.z/cellSz.z));}
	
	
	//cellloc is a particle position - cell idx is going to be floor of each coord
	//just gives value at least and greatest corners of cube
	public myVector getVelAtCell(myVector testLoc) {
		//ctrSzHalfNC == (ctr - (halfNumCell * cellSz))/cellSz
		myVector testLocInFluid = transLocIntoCell(testLoc);
				//new myVector(ctrSzHalfNC,new myVector (testLoc.x/cellSz.x, testLoc.y/cellSz.y, testLoc.z/cellSz.z));		//tmpTestLoc - ctrSzHalfNC;
		
		//double locX = (testLoc(0) - center(0)) / cellSz.x + halfNmCellX,
		//	locY = (testLoc(1) - center(1)) / cellSz.y + halfNmCellY,
		//	locZ = (testLoc(2) - center(2)) / cellSz.z + halfNmCellZ;
		int intLocX = (int)testLocInFluid.x, //idx location in grid
			intLocY = (int)testLocInFluid.y,
			intLocZ = (int)testLocInFluid.z;
		

		double interpX = testLocInFluid.x - intLocX, interpM1X = 1 - interpX,
			   interpY = testLocInFluid.y - intLocY, interpM1Y = 1 - interpY,
			   interpZ = testLocInFluid.z - intLocZ, interpM1Z = 1 - interpZ;
		//bound idx's
		//int tX[2], tY[2], tZ[2];
		//tX[0] = forceIDXBnd(intLocX, s1i[0]),
		//tY[0] = forceIDXBnd(intLocY, s1i[1]),
		//tZ[0] = forceIDXBnd(intLocZ, s1i[2]),
		//tX[1] = forceIDXBnd(intLocX + 1, s1i[0]),
		//tY[1] = forceIDXBnd(intLocY + 1, s1i[1]),
		//tZ[1] = forceIDXBnd(intLocZ + 1, s1i[2]);
	/*
		vector<int> idxs(8); int cnt = 0;
		for (int z = 0; z < 2; ++z) {for (int y = 0; y < 2; ++y) {for (int x = 0; x < 2; ++x) { idxs[cnt++] = IX(tX[x], tY[y], tZ[z]);}}}
	*/	
		int idx000 = IX(forceIDXBnd(intLocX, s1i[0],0),	forceIDXBnd(intLocY, s1i[1],0),	forceIDXBnd(intLocZ, s1i[2],0));
		if(isOOB[idx000]){return new myVector();}
		int idx111 = IX(forceIDXBnd(intLocX + 1, s1i[0],0), forceIDXBnd(intLocY + 1, s1i[1],0),	forceIDXBnd(intLocZ + 1, s1i[2],0));
		if(isOOB[idx111]){return new myVector();}
		//double valx = interpM1X * Vx[idx000] + interpX*Vx[idx111],
		//	   valy = interpM1Y * Vy[idx000] + interpY*Vy[idx111],
		//	   valz = interpM1Z * Vz[idx000] + interpZ*Vz[idx111];

		//int idx = IX(forceIDXBnd((int)((testLoc(0) - center(0)) / cellSz.x + halfNmCellX), s1i[0]),
		//			 forceIDXBnd((int)((testLoc(1) - center(1)) / cellSz.y + halfNmCellY), s1i[1]),
		//			 forceIDXBnd((int)((testLoc(2) - center(2)) / cellSz.z + halfNmCellZ), s1i[2]));

		//return myVector(Vx[idx], Vy[idx], Vz[idx]);
//		myVector res = new myVector(interpM1X * Vx[idx000] + interpX*Vx[idx111], interpM1Y * Vy[idx000] + interpY*Vy[idx111], interpM1Z * Vz[idx000] + interpZ*Vz[idx111]);
//		if(res.magn > 0){
//			pa.pr("idx0 : " + idx000 + " idx1 : " + idx111 + "\tres force : " +  res.toString() + "\ttest loc:"+testLocInFluid.toStrBrf() + "\tpart location normed by cell sz : :" + tmpTestLoc.toStrBrf()+"\n");
//		}
//		return res;
		return new myVector(interpM1X * Vx[idx000] + interpX*Vx[idx111], interpM1Y * Vy[idx000] + interpY*Vy[idx111], interpM1Z * Vz[idx000] + interpZ*Vz[idx111]);
	}//getVelAtCell
	
	
/////////////////////////////////////////////////////////////
///rendering functions
////////////////////////////////////////////////////////////
	private void initGlobeStand(){
		//draw with "shape(globeBase);"
		globeBase = pa.createShape(PConstants.GROUP); 
		globeBase.addChild(buildFloor());
		globeBase.addChild(buildFence(SnowGlobeWin.groundRadius/50.0f));
		globeBase.addChild(buildGlobeBase());	
		
	}//initPieFloor	
	private PShape buildGlobeBase(){
		PShape shRes = pa.createShape(PConstants.GROUP);
		shRes.rotate(PConstants.HALF_PI,-1,0,0);
		myVectorf trans = new myVectorf(0,0,0);
		shRes.addChild(pa.buildCylinder( SnowGlobeWin.groundRadius, gndRad11, gndRadOv6,baseRes,false,false, 3,trans));
		trans._add(new myVectorf(0, gndRadOv6,0)); 
		shRes.addChild(pa.buildCylinder( gndRad11, gndRad11, gndRadOv8,baseRes,false,false,2,trans));	
		trans._add(new myVectorf(0, gndRadOv8,0));
		shRes.addChild(pa.buildCylinder(gndRad12, gndRad12, gndRadOv8,baseRes,true,false,0,trans));
		trans._add(new myVectorf(0, gndRadOv8,0));
		shRes.addChild(pa.buildCylinder(gndRad12,gndRad14, gndRadOv6,baseRes,false,false,3,trans));
		trans._add(new myVectorf(0, gndRadOv6,0));
		shRes.addChild(pa.buildCylinder(gndRad14, gndRad14, gndRadOv6,baseRes,false,false,2,trans));			
	    return shRes;
	}//buildGlobeBase
	
	private PShape buildFloor(){
		PShape shRes = pa.createShape(PConstants.GROUP);
	    int ringIncrVal = 1, oldIncrVal = ringIncrVal;	    
	    float outRad, inRad;
	    //drawPie proto : # of slices around circle, radius, increment amount of drawing for loop, compare value for color change, cosine array, sine array
	    inRad = SnowGlobeWin.groundRadius;
		shRes.rotate(PConstants.HALF_PI,-1,0,0);
	    for(int i=(numGndTarRings-1);i>=0;--i){
	    	int div = 1 << (PApplet.round(PApplet.sqrt(i))+2);
	    	ringIncrVal = pa.numThVals/div;
	    	//System.out.println(" i : " + i + " div : " + div + " ringincr val : " + ringIncrVal + " old val : " +  oldIncrVal + " transition : " + (ringIncrVal!=oldIncrVal));
	    	outRad = inRad;//decreasing values
	    	inRad = i * radRat;	    	
			shRes.addChild(buildPieRing(outRad,  inRad, ringIncrVal, 0, ringIncrVal!=oldIncrVal, i%2==0));
			oldIncrVal = ringIncrVal;
	    }
		return shRes;
	}	
	
	/**
	*  builds the concentric rings of the "chessboard" into a pshape
	*  @param pieRadOutside the radius of this piePiece on the outside
	*  @param pieRadInside the radius of this piePiece on the inside
	*  @param incrLoop the amount to increment the loop each time - increment by more as we get closer to the center
	*  @param y height
	*  @param draw5th whether to draw a 5th vertex to align transition rings from lower to higher count of tiles
	*  @paramm whether to start with a black or a white tile - always want to alternate
	*/	
	private PShape buildPieRing(float pieRadOutside, float pieRadInside, int incrLoop, float y, boolean draw5th, boolean initDrawBlack){
		PShape shRes = pa.createShape(PConstants.GROUP),sh;
	    int ii;    
	    int incr2  = (int)(incrLoop/2);
	    boolean drawBlack = initDrawBlack;
	    for (int i = 0; i < pa.numThVals; i += incrLoop) {
	    	ii = (i+incrLoop) % pa.numThVals;
	    	sh = pa.createShape();
	    	sh.beginShape();
	    	if (drawBlack){ 		setClrShinyBlack(sh);   } 
	    	else{					setClrFlatWhite(sh);    }
	    	drawBlack = !drawBlack;
	    	sh.stroke(128,100,0,255);
	    	sh.strokeWeight(1.0f);
	    	sh.vertex ((pa.cosAra[ii] * pieRadOutside), y,	(pa.sinAra[ii] * pieRadOutside));//anchor point for vertex
        	if (draw5th){  
        		int idx2 = (i+incr2) % pa.numThVals;
        		sh.vertex ((pa.cosAra[idx2] * pieRadOutside), y, (pa.sinAra[idx2] * pieRadOutside));      }
        	sh.vertex ((pa.cosAra[i] * pieRadOutside),y,(pa.sinAra[i] * pieRadOutside));
        	sh.vertex ((pa.cosAra[i] * pieRadInside),y,(pa.sinAra[i] * pieRadInside));
        	sh.vertex ((pa.cosAra[ii] * pieRadInside),y,(pa.sinAra[ii] * pieRadInside));
        	//sh.noStroke();
        	sh.endShape(pa.CLOSE);
			shRes.addChild(sh);
	    }//for
	    return shRes;
	}//buildPieRing method	
	
	private PShape buildFence(float cylSmHght){
		float smHght = cylSmHght;
		PShape shRes = pa.createShape(PConstants.GROUP);
		int ii, incr = pa.numThVals/baseRes;
		myVectorf[] verts = new myVectorf[]{new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0)};
		shRes.rotate(PConstants.HALF_PI,1,0,0);
		for(int j=0;j<numFencePosts;++j){
	    	if(j%3 == 0) { shRes.addChild(pa.buildCylinder( SnowGlobeWin.groundRadius, SnowGlobeWin.groundRadius, cylSmHght, baseRes, false,false,3,new myVectorf(0,smHght-cylSmHght,0)));} 
		    else{		    	
				for (int i = 0; i < pa.numThVals; i+=incr) {
		    		if ((i/incr) % (numFencePosts*postWidth) < (postWidth+1)){
		    			ii = (i+incr) % pa.numThVals;
						//need to send precalced normal and vertex arrays
						//idx 0,3 are normal values - use half height diff for titled cylinders
		    			verts[0].set((pa.cosAra[i] * SnowGlobeWin.groundRadius), 0, (pa.sinAra[i] * SnowGlobeWin.groundRadius));
		    			verts[1].set((pa.cosAra[i] * SnowGlobeWin.groundRadius), smHght, (pa.sinAra[i] * SnowGlobeWin.groundRadius));
		    			verts[2].set((pa.cosAra[ii] * SnowGlobeWin.groundRadius), smHght, (pa.sinAra[ii] * SnowGlobeWin.groundRadius));
		    			verts[3].set((pa.cosAra[ii] * SnowGlobeWin.groundRadius), 0, (pa.sinAra[ii] * SnowGlobeWin.groundRadius));
		    			shRes.addChild(pa.buildVShape(verts, 3));
		    		}
				}		    	
		    }
	    	smHght += cylSmHght;	    	  
	    }
	    return shRes;
	}//drawFence
	
	public void drawGlobe(){
		pa.pushMatrix();pa.pushStyle();
		    pa.setColorSnowGlobeBall(); 
		    pa.sphereDetail(120);
		    pa.translate(center.x,center.y,center.z);
		    pa.unSetCamOrient();
		    pa.rotate(PConstants.HALF_PI, 1,0,0);
		    pa.sphere((float)(SnowGlobeWin.groundRadius*1.08));
		    pa.setColorClearWhite(); 
		    pa.sphere((float)(SnowGlobeWin.groundRadius*1.09));
	    pa.popStyle();pa.popMatrix();
	}
	
	private void setClrShinyBlack(PShape sh){		setShColorVals(sh, 0,0,0, 255, 255, 255, 255,300.0f);}
	private void setClrFlatWhite(PShape sh){		setShColorVals(sh,255,255,255, 255,55, 55, 55, 1.0f);}

	private void setShColorVals(PShape sh, int f_r, int f_g, int f_b, int f_a, 
		  					int sp_r, int sp_g, int sp_b, 
		  					float si){		  
	    sh.fill(f_r, f_g, f_b, f_a);
	    sh.specular(sp_r, sp_g, sp_b);
	    sh.shininess(si);
	}//setColorVals
	
}//snoGlobe class

class myFluidBndObj {
	public final int idx;//, x, y, z;
	public myVectorf norm, loc, frcAdded;
	public double mag, oneMMag;					//0-1, pct of cell within boundary

	//public myFluidBndObj(int _x, int _y, int _z, int _idx, myVectorf _normToCtr, double _mag) {
	public myFluidBndObj(int _idx, myVectorf _normToCtr, myVectorf _loc,  double _mag) {
		//x=_x;y=_y;z=_z;
		idx=_idx;
		norm = new myVectorf(_normToCtr);frcAdded = new myVectorf();
		loc = new myVectorf(_loc);
		mag = _mag;
		oneMMag = 1 - mag;
	}
	
}//myFluidBndObj class


