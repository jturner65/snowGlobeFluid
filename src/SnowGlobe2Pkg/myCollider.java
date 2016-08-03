package SnowGlobe2Pkg;

import SnowGlobe2Pkg.SnowGlobeWin.CollisionType;
import SnowGlobe2Pkg.SnowGlobeWin.SolverType;


public abstract class myCollider {
	protected static SnowGlobeWin pa;
	protected static mySnowGlobeWin win;
	public static int ID_gen = 0;
	public int ID;
	public String name;
	public SnowGlobeWin.CollisionType colType;
	
	public myVectorf drawLoc;			//drawn location of this collider - only different from center if needed for display purposes

	public double Krest,			//coefficent of restitution - how much bounce do we have :1 total bounce, 0 no bounce.  multiply this against part's Velperp
				muFrict;         //friction coefficient

	public static final int
		NoCol = 0,					//0 if no collision, 
		BrchCol = 1,				//1 if collision via breach - push back to legal position, address vel and frc concerns
		NextCol = 2,				//2 if collision next timestep  
		CntctCol = 3;				//3 if need to counter force due to contact - within some epsilon of part radius distance from collider surface

		
	
	public static final double partRad = .04;
	
	public myCollider(SnowGlobeWin _p, mySnowGlobeWin _win, String _n, myVectorf _drawLoc, SnowGlobeWin.CollisionType _colType) {
 		pa = _p;
  		win = _win;
  		ID = ID_gen++;
  		name = new String(_n);
  		drawLoc = _drawLoc;
  		colType = _colType;
	}
	
	//checks if particle location + deltaT partVel will cause a collision, or if particle is in contact without
	//a normal-dir component of velocity, in which case it will modify the force
	//0 vec if no collision, 1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact	
	public abstract int checkCollision(double deltaT, myParticle part);

	//if particle has breached collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	public abstract void handleCollision(myParticle part, int res);
	
}//myCollider

class cylinderCollider extends myCollider{
	public myVector center,			//actual location of this collider in the world
		dims,						//x,z : radii along ciruclar axes, around orient axis. y : height along orient axis above center
		orient;						//should be normalized vector pointing in direction of cylinder axis		
	public double[] minMaxDims;		//minimum dist from center to still be inside (in case of ellipse, to minimize calculations) min idx 0, max idx 1

	public boolean  intRefl;			//internal reflections? for sphere or cylinder (collide on inside instead of outside)
	
	public cylinderCollider(SnowGlobeWin _p, mySnowGlobeWin _win, String _n, myVectorf _drawLoc, myVector _ctr, myVector _dims, myVector _orient,  boolean _intRefl) {
		super(_p, _win, _n, _drawLoc, CollisionType.CYLINDER);
		intRefl = _intRefl;
		center = _ctr; dims = _dims; orient = _orient;
		findMinMaxDims();
	}

	//finds minimum and maximum value of radius for cylinder, to speed up calculations of collisions
	public void findMinMaxDims() {//idx 0 is min, idx 1 is max
		minMaxDims = new double[4];
		minMaxDims[0] = (dims.z < dims.x ? dims.z : dims.x);
		minMaxDims[1] = (dims.z < dims.x ? dims.x : dims.z);
		minMaxDims[2] = minMaxDims[0] * minMaxDims[0];	//sq min rad
		minMaxDims[3] = minMaxDims[1] * minMaxDims[1];	//sq max rad
	}	

	@Override
	//0 vec if no collision, 1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	public int checkCollision(double deltaT, myParticle part) {
		myVector partPos = part.aPosition[part.curIDX];
		//TODO need to perform quick calc
		
		double hfDelT2 = .5 * deltaT * deltaT;
		myVector multPartFAcc = myVector._mult(part.aForceAcc[part.curIDX], hfDelT2),		
		//forward integrate for estimate of future
				partVelPoint = myVector._add(myVector._mult(part.aVelocity[part.curIDX],deltaT), multPartFAcc),
				partMovePoint = myVector._add(partPos, partVelPoint);	//potential movement point for next turn of movement, to see if next turn of movement will hit wall
		//TODO fix this to handle cylinder
		myVector partMvCtr = new myVector(center,partMovePoint), partLocCtr = new myVector(center,partPos);
			if (((partMvCtr.magn < (minMaxDims[0] + partRad)) && (!intRefl)) ||					//current location is breach 
			((partMvCtr.magn > (minMaxDims[1] - partRad)) && (intRefl))) {
				if (((partLocCtr.magn < (minMaxDims[0] + partRad)) && (!intRefl)) ||					//current location is breach 
					((partLocCtr.magn > (minMaxDims[1] - partRad)) && (intRefl))) {
					return BrchCol;
				}
				else {
					return NextCol;
				}
			}
		//find point on surface of sphere inline with center and partlocation
		myVector sNormPartP = getCylinderNormal(partPos);				//normal through current point and center, in direction of collision surface
		myVector partSpherePnt =  myVector._add(center, myVector._mult(sNormPartP,-snoGlobe.snowGlobRad));			//point on ellipsoid surface colinear with center and particle move point
		double dist2wall = myVector._sub(partSpherePnt, partPos).magn, distFromWallChk = dist2wall - partVelPoint.magn;
		if (distFromWallChk > pa.epsVal) { return NoCol; }
		else if (distFromWallChk > -pa.epsVal) { return CntctCol; }
		else { return BrchCol; }

	}//checkCollision	
	
	public myVector getCylinderNormal(myVector _loc) {//get normal at a particular location - no matter where inside or outside of sphere, normal built from this point and center will point in appropriate dir
		//if sphere, normal will be either pointing out or in, colinear with line from center to _loc
		myVector normDir = myVector._sub(center,_loc);//either point into center if internal reflections or point out of center if not
		double mult = ((intRefl) ? 1 : -1);
		normDir._mult(mult);
		normDir._normalize();
		return normDir;
	}//getNormal

	@Override
	//if particle has breached collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	public void handleCollision(myParticle part, int res) {
		//myVector partPos = part.position.peekFirst();
		myVector partPos = part.aPosition[part.curIDX];
		//myVector partVel = part.velocity.peekFirst();
		myVector partVel = part.aVelocity[part.curIDX];
		myVector cylNormal = getCylinderNormal(partPos);

		if (res == 2) {//if close to intersection with sphere boundary
			myVector[] partVelComp = pa.getVecFrame(partVel, cylNormal);
			//partVelComp[0] *= (-1 * Krest);//reverse direction of normal velocity
			partVel.set(myVector._add(myVector._mult(partVelComp[0],(-1 * Krest)), partVelComp[1]));//should change dir of velocity, decrease tangent velocity for friction
		}//if about to hit collider

		else if (res == 1) {//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
			double distFromBreach = myVector._sub(partPos, center).magn - (snoGlobe.snowGlobRad - partRad);
			//cout<<"dist from breach "<<distFromBreach<<endl;
			if (((intRefl) && ((distFromBreach) < 0)) || ((!intRefl) && ((distFromBreach) > 0))) {}//cout<<"breach error, not on wrong side of sphere"<<endl;}
			else {//forcibly move particle to just a bit on the right side of the collider, reverse velocity
				distFromBreach *= (1.1);//move slightly more than breach amount
				myVector newPos = myVector._add(partPos,myVector._mult(cylNormal,distFromBreach));//move back into sphere
				partPos.set(newPos);
				//part.position.peekFirst().set(newPos);
				myVector[] partVelComp = pa.getVecFrame(partVel, cylNormal);
				//partVelComp[0] *= -1;//reverse direction of normal velocity
				partVel.set(myVector._mult(myVector._add(partVelComp[0],partVelComp[1]),-1));//should change dir of velocity, for sphere zeroing tangent velocity
			}
		}//if 1
		else if (res == 3) //|| (res == 2) || (res == 1)) 
		{//tangent, get forceAcc and add -(forcecomponent in normal dir)
			myVector[] partAccComp = pa.getVecFrame(part.aForceAcc[part.curIDX], cylNormal);
			partAccComp[0]._mult( -1 * Krest);//reverse direction of normal accel
			part.applyForce(myVector._add(partAccComp[0], partAccComp[1]));
		}//tangent
	}//handlePlanarBreach	
	
}//cylinderCollider


class sphereCollider extends myCollider{
	public myVector center,			//actual location of this collider, w/respect to particles - used for sphere
		radius;						//radius around center for ellipsoid (in each direction for x,y,z)
	public double[] minMaxRadius;		//minimum dist from center to still be inside (in case of ellipse, to minimize calculations) min idx 0, max idx 1


	public boolean  intRefl;			//internal reflections? for sphere (collide on inside)
//	myCollider(string _n, const Eigen::Vector3d& _dr, const Eigen::Vector3d& _ctr, const Eigen::Vector3d& _rad, bool _inRefl) :							//sphere collider
//		ID(++ID_gen), name(_n), colType(SPHERE), drawLoc(_dr), center(_ctr), radius(_rad), minMaxRadius(2), planeNormal(), verts(), peq(4), intRefl(_inRefl), Krest(1) {
//		initCollider();
//	}

	public sphereCollider(SnowGlobeWin _p, mySnowGlobeWin _win, String _n, myVectorf _drawLoc, myVector _ctr, myVector _rad, boolean _intRefl) {
		super(_p, _win, _n, _drawLoc, CollisionType.SPHERE);
		intRefl = _intRefl;
		center = _ctr; radius = _rad;
		findMinMaxRadius();
	}

	//finds minimum and maximum value of radius for ellipsoid sphere, to speed up calculations of collisions
	public void findMinMaxRadius() {
		minMaxRadius = new double[5];
		minMaxRadius[0] = pa.min3(radius.z, radius.x, radius.y);
		minMaxRadius[1] = pa.max3(radius.z, radius.x, radius.y);
		minMaxRadius[2] = minMaxRadius[0] * minMaxRadius[0];	//sq min rad
		minMaxRadius[3] = minMaxRadius[1] * minMaxRadius[1];	//sq max rad
		minMaxRadius[4] = minMaxRadius[0] - pa.tenSFRads;		//min dist to ignore particle collision - min radius - 10x particle radius
	}	

	@Override
	//0 vec if no collision, 1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	public int checkCollision(double deltaT, myParticle part) {
		myVector partLocVec = new myVector(center,part.aPosition[part.curIDX]), vecToCtr = myVector._normalize(partLocVec);			//vector from center to particle position to get dist sphere wall
		//far from plane - no need to check further in collision detection
		double distFromCtr = partLocVec.magn + pa.snowFlakeRad;
		if (distFromCtr < minMaxRadius[4]){ return NoCol;}						// more than 10x part radii from sphere surface - no collision assumed possible
		double distFromCtrDiff  = distFromCtr - minMaxRadius[1];			//compare to largest dimension of ellipsoid - if still positive, then definite breach
		if (distFromCtrDiff > pa.epsVal) { return BrchCol; }						//immediate collision - breached plane by more than eps + snowflakerad

		double spdInRadNormDir = part.aVelocity[part.curIDX]._dot(vecToCtr),				//velocity in direction of particle-from-center vector
				accInRadNormDir = part.aForceAcc[part.curIDX]._dot(vecToCtr)/part.mass;		//acc in direction of particle-from-center vector
	
		if((spdInRadNormDir > 0) && (accInRadNormDir > 0)){ return NoCol;}						//not moving toward plane or accelerating toward plane, so no collision possible this time step
		//by here, within col dist of plane, and moving toward, or tangent to, plane - predict motion
		double velPartInRadNormDir = spdInRadNormDir * deltaT,
				accPartInRadNormDir = .5 * deltaT * deltaT * accInRadNormDir,
				distWillMove = velPartInRadNormDir + accPartInRadNormDir;
		
		myVector v2ctrORad = myVector._elemDiv(vecToCtr, radius);
		
		double a = v2ctrORad._dot(v2ctrORad), b = v2ctrORad._dot(center), c = center._dot(center) -1, ta = 2*a, discr1 = Math.pow(((b*b) - (2*ta*c)),.5), 
				t1 = (-1*b + discr1)/(ta), t2 = (-1*b - discr1)/(ta);
		//set the t value of the intersection to be the minimum of these two values (which would be the edge closest to the eye/origin of the ray)
		double t = Math.max(t1,t2);
		
		if(distFromCtr + distWillMove > t){//will move further from center than location of intersection of vector from center
			return NextCol;
		}
		return NoCol;

//		//calc t so that normalized partLocVec collides with sphere/ellipsoid wall.  if t > len(partLocVec) then no collision
//		
//		
//		
//		
//		
//		myVector partPos = part.aPosition[part.curIDX],				
//				partLocCtr = new myVector(center,partPos);
//		//check if partloc is very far from wall
//		if(((partLocCtr.magn < .95f * minMaxRadius[0]) && intRefl) //less than 90% of the sphere's minimum radius or 111% of max radius
//			|| ((partLocCtr.magn > 1.056f * minMaxRadius[1]) && !intRefl)){return NoCol;}
//		
//		double hfDelT2 = .5 * deltaT * deltaT;
//		myVector multPartFAcc = myVector._mult(part.aForceAcc[part.curIDX], hfDelT2),
//				partVelPoint = myVector._add(myVector._mult(part.aVelocity[part.curIDX],deltaT), multPartFAcc),
//				partMovePoint = myVector._add(partPos, partVelPoint);	//potential movement point for next turn of movement, to see if next turn of movement will hit wall
//
//		myVector partMvCtr = new myVector(center,partMovePoint);
//			if (((partMvCtr.magn < (minMaxRadius[0] + partRad)) && (!intRefl)) ||					//current location is breach 
//			((partMvCtr.magn > (minMaxRadius[1] - partRad)) && (intRefl))) {
//				if (((partLocCtr.magn < (minMaxRadius[0] + partRad)) && (!intRefl)) ||					//current location is breach 
//					((partLocCtr.magn > (minMaxRadius[1] - partRad)) && (intRefl))) {
//					return BrchCol;
//				}
//				else {
//					return NextCol;
//				}
//			}
//		//find point on surface of sphere inline with center and partlocation
//		myVector sNormPartP = getSphereNormal(partPos);				//normal through current point and center, in direction of collision surface
//		myVector partSpherePnt =  myVector._add(center, myVector._mult(sNormPartP,-snoGlobe.snowGlobRad));			//point on ellipsoid surface colinear with center and particle move point
//		double dist2wall = myVector._sub(partSpherePnt, partPos).magn, distFromWallChk = dist2wall - partVelPoint.magn;
//		if (distFromWallChk > pa.epsVal) { return NoCol; }
//		else if (distFromWallChk > -pa.epsVal) { return CntctCol; }
//		else { return BrchCol; }
	}//checkCollision		
	
	public myVector getSphereNormal(myVector _loc) {//get normal at a particular location - no matter where inside or outside of sphere, normal built from this point and center will point in appropriate dir
		//if sphere, normal will be either pointing out or in, colinear with line from center to _loc
		myVector normDir = myVector._sub(center,_loc);//either point into center if internal reflections or point out of center if not
		double mult = ((intRefl) ? 1 : -1);
		normDir._mult(mult);
		normDir._normalize();
		return normDir;
	}//getNormal
	
	@Override
	//if particle has breached planar collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	public void handleCollision(myParticle part, int res) {
		myVector partPos = part.aPosition[part.curIDX], partVel = part.aVelocity[part.curIDX], sphereNormal = getSphereNormal(partPos);

		if (res == 2) {//if close to intersection with sphere boundary
			myVector[] partVelComp = pa.getVecFrame(partVel, sphereNormal);
			//partVelComp[0] *= (-1 * Krest);//reverse direction of normal velocity
			partVel.set(myVector._add(myVector._mult(partVelComp[0],(-1 * Krest)), partVelComp[1]));//should change dir of velocity, decrease tangent velocity for friction
		}//if about to hit collider

		else if (res == 1) {//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
			double distFromBreach = myVector._sub(partPos, center).magn - (snoGlobe.snowGlobRad - partRad);
			//cout<<"dist from breach "<<distFromBreach<<endl;
			if (((intRefl) && ((distFromBreach) < 0)) || ((!intRefl) && ((distFromBreach) > 0))) {}//cout<<"breach error, not on wrong side of sphere"<<endl;}
			else {//forcibly move particle to just a bit on the right side of the collider, reverse velocity
				distFromBreach *= (1.1);//move slightly more than breach amount
				myVector newPos = myVector._add(partPos,myVector._mult(sphereNormal,distFromBreach));//move back into sphere
				partPos.set(newPos);
				myVector[] partVelComp = pa.getVecFrame(partVel, sphereNormal);
				//partVelComp[0] *= -1;//reverse direction of normal velocity
				partVel.set(myVector._mult(myVector._add(partVelComp[0],partVelComp[1]),-1));//should change dir of velocity, for sphere zeroing tangent velocity
			}
		}//if 1
		else if (res == 3) //diminish all force and velocity in normal dir 
		{//tangent, get forceAcc and add -(forcecomponent in normal dir)
			myVector[] partAccComp = pa.getVecFrame(part.aForceAcc[part.curIDX], sphereNormal);
			partAccComp[0]._mult( -1 * Krest);//reverse direction of normal accel
			part.applyForce(myVector._add(partAccComp[0], partAccComp[1]));
		}//tangent
	}//handlePlanarBreach	

}//sphereCollider

class planeCollider extends myCollider{
	//plane
	public myVector planeNormal;		//normal of this collider, if flat plane
	public myVector[] verts;	//vertices of this object, if flat plane
	public double[] peq;		//plane equation values
	
	public planeCollider(SnowGlobeWin _p, mySnowGlobeWin _win, String _n, myVectorf _drawLoc, myVector[] _verts) {
		super(_p, _win, _n, _drawLoc, CollisionType.FLAT);
		verts = _verts;
		buildPlaneNorm();
		findPlaneEQ();
	}
	
	//determines the equation coefficients for a planar collider
	public void findPlaneEQ() {
		//Ax + By + Cz + D = 0
		peq = new double[4];
		peq[0] = planeNormal.x;		//A == norm.X
		peq[1] = planeNormal.y;		//B == norm.y
		peq[2] = planeNormal.z;		//C == norm.z
		peq[3] = -planeNormal._dot(verts[0]);//- ((peq[0] * verts[0].x) + (peq[1] * verts[0].y) + (peq[2] * verts[0].z));		//D
	}

	//build normal of planar object
	public void buildPlaneNorm() {
		myVector P0P1 = myVector._sub(verts[1], verts[0]);
		myVector P1P2 = myVector._sub(verts[1], verts[2]);
		//planeNormal = P0P1._cross(P1P2);
		planeNormal = P1P2._cross(P0P1);
		planeNormal._normalized();
	}//buildNorm

	@Override
	//0 if no collision, 
	//1 if collision via breach - push back to legal position, address vel and frc concerns
	//2 if collision next timestep  
	//3 if need to counter force due to contact - within some epsilon of part radius distance from collider surface
	public int checkCollision(double deltaT, myParticle part) {
		myVector partLocVec = new myVector(verts[0],part.aPosition[part.curIDX]);			//vector from point on plane to particle position to get dist from plane
		//far from plane - no need to check further in collision detection
		double distFromPlane = partLocVec._dot(planeNormal) - pa.snowFlakeRad; 				//distance edge of particle is from plane
		if (distFromPlane  > pa.tenSFRads){ return NoCol;}									//if further away from plane than 10 * particle radius then no collision this or next cycle
		if (distFromPlane  < -pa.epsVal) { return BrchCol; }						      	//immediate collision - breached plane by more than eps*rad
		//dist between epsVal and 10 snoflake rads - possible collision next cycle
		double spdInPlaneNormDir = part.aVelocity[part.curIDX]._dot(planeNormal),					//velocity in direction of plane normal - speed toward plane is negative of this
				accInPlaneNormDir = part.aForceAcc[part.curIDX]._dot(planeNormal)/part.mass;		//acc in dir of plane normal - acc toward plane is negative of this
		if((spdInPlaneNormDir > 0) && (accInPlaneNormDir > 0)){ return NoCol;}						//not touching plane, moving toward plane or accelerating toward plane, so no collision possible this time step
		if (distFromPlane  < pa.epsVal) { return CntctCol; }								//contact - address forces
		//by here, within col dist of plane, and moving toward, or tangent to, plane - predict motion
		double hfDelT2 = .5 * deltaT * deltaT,
				velPartInPlaneDir = spdInPlaneNormDir * deltaT,
				accPartInPlaneDir = .5 * deltaT * deltaT * accInPlaneNormDir,
				distWillMove = velPartInPlaneDir + accPartInPlaneDir;		
		if(distFromPlane < distWillMove){//particle is closer to plane than how far it is going to move next cycle - will breach after integration
			return NextCol;
		}
		return NoCol;

	}//checkCollision	
	
	//if particle has breached planar collider somehow, move particle along appropriate normal direction until not breached anymore, and check particle's velocity(and reverse if necessary)
	//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
	@Override
	public void handleCollision(myParticle part, int res) {
		myVector partPos = part.aPosition[part.curIDX];
		double distFromBreach = myVector._sub(partPos, verts[0])._dot(planeNormal);
		myVector partVel = part.aVelocity[part.curIDX], partFrc = part.aForceAcc[part.curIDX];		
		
		//1 if collision via breach, 2 if collision next timestep, 3 if need to counter force due to contact
		if ((res == 2) || (partVel._dot(planeNormal) < 0)) {//going to breach next time step, or have breached and velocity is still going down - swap velocity direction
			myVector[] partVelComp = pa.getVecFrame(partVel, planeNormal);
			//partVelComp[0] *= (-1 * Krest);//reverse direction of normal velocity
			partVel.set(  myVector._add(myVector._mult(partVelComp[0],(-1 * Krest)), partVelComp[1]));//should change dir of velocity
			// part.color = myVector(1,0,0);
			if (part.solveType == SolverType.VERLET) { handleVerletCol(part); }//handle reflection/velocity change by swapping old and new positions - need to scale by krest
		}//if breached and velocity going away from wall

		else if (res == 1) {           //immediate breach, swap position
			//p.color = myVector(0, 0, 0);
			if (distFromBreach > 0) {}//cout<<"breach error, not on wrong side of plane"<<endl;}
			else if (part.solveType == SolverType.VERLET) { handleVerletCol(part); }	//handle reflection/velocity change by swapping old and new positions - need to scale by krest			
			else {//forcibly move particle to just a bit on the right side of the collider
				distFromBreach *= -(2.001);
				myVector newPos = myVector._add(partPos, myVector._mult(planeNormal,(distFromBreach + pa.epsVal)));   //reflect position up from plane by slightly more than breach amount
				//if(p.getSolveType() == GROUND){cout<<"dist from breach "<<distFromBreach<<" old position: "<<partPos<<" new position : "<<newPos<<endl;}
				//part.position.peekFirst().set( newPos);
				partPos.set( newPos);
				myVector[] partAccComp = pa.getVecFrame(partFrc, planeNormal);
				myVector frcTanDir = myVector._normalize(partAccComp[1]),
						//TODO fix this stuff - friction is not working correctly
						tanForce = myVector._mult(frcTanDir,-muFrict * (partAccComp[0]._dot(planeNormal)))
				;
				partFrc.set(partAccComp[0]) ;
				partVel.set(0,0,0);// = myVector(0,0,0);//partVelComp[0];//+partVelComp[1];//should change dir of velocity
			}
		}//if 1

		else if (res == 3) {          //contact
			if (part.solveType == SolverType.VERLET) { handleVerletCol(part); }//handle reflection/velocity change by swapping old and new positions - need to scale by krest
		}

		if ((res == 3) || (res == 2) || (res == 1)) {                             //any contact criteria - swap normal force direction
			myVector[] partAccComp = pa.getVecFrame(partFrc, planeNormal);
			partAccComp[0]._mult(-1);//reverse direction of normal acc
			//part.forceAcc.peekFirst().set(myVector._add(partAccComp[0], partAccComp[1]));
			partFrc.set(myVector._add(partAccComp[0], partAccComp[1]));
		}//tangent
	}//handlePlanarBreach

	public void handleVerletCol(myParticle p) {
		myVector tmpOldPos = p.aPosition[p.curIDX], 
				tmpNewPos = p.aOldPos[p.curIDX];         //swapped already
		myVector colPt = myVector._mult(myVector._add(tmpOldPos, tmpNewPos), .5);
		double krTmp = ((1 - Krest) * .25) + Krest;
		myVector tmpOldSub = myVector._sub(tmpOldPos, colPt), colNDotVec = myVector._mult(planeNormal, colPt._dot(planeNormal)),
				tmpNewSub = myVector._sub(tmpNewPos, colPt);
		
		tmpOldPos = myVector._mult(myVector._add(myVector._sub(myVector._mult(myVector._sub(tmpOldPos,colNDotVec),2), tmpOldSub), colPt),krTmp);
		tmpNewPos = myVector._add(myVector._sub(myVector._mult(myVector._sub(tmpOldPos,colNDotVec),2), tmpNewSub),colPt);
		p.aPosition[p.curIDX].set(tmpNewPos);
		p.aOldPos[p.curIDX].set(tmpOldPos);
	}
}//planeCollider


//an object to hold the information about a collision 
class collisionEvent{
	
	
	
}

//
//class boxCollider extends myCollider{
//	//check if inside box of planes bounded by 8 verts - treat like 6 bounding planes 
//	public planeCollider[] box;
//
//	public boxCollider(SnowGlobeWin _p, mySnowGlobeWin _win, String _n, myVectorf _drawLoc, myVector[] _verts) {
//		super(_p, _win, _n, _drawLoc, CollisionType.BOX);
//		box = new planeCollider[6];
//		for(int i =0;  i<3;++i){//0->3, 1->4, 2->5 are parallel planes, idxs 0-3 are one plane, 4-7 are parallel plane
//			//String _n, myVectorf _drawLoc, myVector[] _verts)
//			box[i] = new planeCollider(_p, _win, _n+"Side_"+i+"a",);
//		}
//		// TODO Auto-generated constructor stub
//	}
//
//	@Override
//	public int checkCollision(double deltaT, myParticle part) {
//		// TODO Auto-generated method stub
//		return NoCol;
//	}
//
//	@Override
//	public void handleCollision(myParticle part, int res) {
//		// TODO Auto-generated method stub
//		
//	}
//	
//}//boxCollider