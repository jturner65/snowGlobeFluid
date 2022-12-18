package SnowGlobe2Pkg;

public abstract class myForce {
	protected static SnowGlobeWin pa;
	public static int ID_gen = 0;
	public int ID;
	public String name;
	public double constVal;				//multiplicative constant to be applied to mass to find force
	public double constVal2;
	public myVector constVec;				//vector constant quantity, for use with gravity
	public SnowGlobeWin.ForceType ftype;

	public myForce(SnowGlobeWin _p, String _n, double _k1, double _k2, myVector _constVec, SnowGlobeWin.ForceType _t){
  		pa = _p; 		ID = ++ID_gen;
		name = new String(_n);
		constVal = _k1; 
		constVal2 = _k2;
		constVec = _constVec;		//torque-result force
		ftype = _t;
	}
	public myForce(SnowGlobeWin _p,String _n, double _k1, double _k2) {this(_p,_n, _k1, _k2, new myVector(), SnowGlobeWin.ForceType.DAMPSPRING);}
	public myForce(SnowGlobeWin _p,String _n, double _k) {this(_p,_n, _k * (_k>0 ? 1 : -1), 0, new myVector(), (_k>0) ? SnowGlobeWin.ForceType.REPL : SnowGlobeWin.ForceType.ATTR); ID = -1;}
	
	public abstract myVector[] calcForceOnParticle(myParticle _p1, myParticle _p2, double d);// {S_SCALAR,S_VECTOR, ATTR, SPRING};
	@Override
	public String toString(){return "Force Name : " + name + " ID : " + ID + " Type : " + pa.ForceType2str[ftype.getVal()];}
}//myForce class

class mySclrForce extends myForce{
	// "scalar" force here means we derive the force by a particle-dependent scalar value, in this case mass against gravity vec 
	public mySclrForce(SnowGlobeWin _p,String _n, myVector _G) { super(_p,_n, 0 ,0, new myVector(_G), SnowGlobeWin.ForceType.S_SCALAR);}	//	

	@Override
	//array returns up to 2 forces, one on p1, one on p2
	public myVector[] calcForceOnParticle(myParticle _p1, myParticle _p2, double d) {
		myVector[] result = new myVector[]{new myVector(),new myVector()};
		result[0] = myVector._mult(constVec,_p1.mass);
		return result;
	}
	@Override
	public String toString(){return super.toString() + "\tForce Vector :  " + constVec.toString();}
	
}//mySclrForce - scalar body-specific multiple of vector force

class myVecForce extends myForce{
	//vector here means we derive the force as a particle-dependent vector value, like velocity, against some scalar kd
	public myVecForce(SnowGlobeWin _p,String _n, double _k) { super(_p,_n,_k,0, new myVector(), SnowGlobeWin.ForceType.S_VECTOR);}		//if drag, needs to be negative constant value	

	@Override
	public myVector[] calcForceOnParticle(myParticle _p1, myParticle _p2, double d) {
		myVector[] result = new myVector[]{new myVector(),new myVector()};
		result[0] = myVector._mult(_p1.aVelocity[_p1.curIDX], constVal);//vector here means we derive the force as a particle-dependent vector value, velocity, against some scalar kd 
		return result;
	}
	@Override
	public String toString(){return super.toString() + "\tForce Scaling Constant :  " + String.format("%.4f",constVal);}
	
}//myVecForce - vector body-specific quantity multiplied by scalar constant

class my2bdyForce extends myForce{
	//attractive/repulsive force
	public my2bdyForce(SnowGlobeWin _p,  String _n, double _k,  SnowGlobeWin.ForceType _t) {
		super(_p, _n, _k, 0, new myVector(), _t);
	}
	public my2bdyForce(SnowGlobeWin _p,  String _n, double _k) {//passed k > 0 is repulsive force, k < 0 is attractive force
		this(_p, _n, Math.abs(_k), (_k>0) ? SnowGlobeWin.ForceType.REPL : SnowGlobeWin.ForceType.ATTR);
	}
	@Override
	public myVector[] calcForceOnParticle(myParticle _p1, myParticle _p2, double d) {
		myVector[] result = new myVector[]{new myVector(),new myVector()};
		myVector vecL;
		vecL = new myVector(_p2.aPosition[_p2.curIDX],_p1.aPosition[_p1.curIDX]);//vector from 2 to 1
		if (vecL.magn > pa.epsVal) {		
			double m1 = _p1.mass, m2 = _p2.mass;
			myVector lnorm = myVector._normalize(vecL);			//unitlength vector of l
			double fp = constVal * m1 * m2 / (vecL.sqMagn);		//from 2 to 1 if constVal > 0 (repulsive force)
			result[0] = myVector._mult(lnorm, fp);				//force applied to p1
			result[1] = myVector._mult(lnorm, -fp);				//force applied to p2
		}//only add force if magnitude of distance vector is not 0
		return result;
	}	
	@Override
	public String toString(){return super.toString() + "\tForce Scaling Constant :  " + String.format("%.4f",constVal);}	
}

//myVector[] result = new myVector[]{new myVector(),new myVector()};
//myVector vecL, v_l;
//switch (force.ftype) {
//	case S_VECTOR:	{
//		result[0] = myVector._mult(_p1.velocity.peekFirst(), force.constVal);//vector here means we derive the force as a particle-dependent vector value, velocity, against some scalar kd 
//		break; }
//	case ATTR: {//attractor, uses two particles, 1st constant, 
//		vecL = new myVector(_p2.position.peekFirst(),_p1.position.peekFirst());//from 2 to 1
//		if (vecL.magn > pa.epsVal) {		//attractor force - constVal (negative) * m1 * m2 * lnorm/ lmag*lmag
//			double m1 = _p1.mass, m2 = _p2.mass;
//			myVector lnorm = myVector._normalize(vecL);//unitlength vector of l
//			double fp = -1 * force.constVal * m1 * m2 / (vecL.sqMagn);
//			result[0] = myVector._mult(lnorm,fp);			
//			result[1] = myVector._mult(lnorm, -fp);
//		}//only add force if magnitude of distance vector is not 0
//		break; }
//	case REPL: {//repulsive force, uses two particles, 1st constant, opposite sign as attractor 
//		vecL = new myVector(_p2.position.peekFirst(),_p1.position.peekFirst());//from 2 to 1
//		if (vecL.magn > pa.epsVal) {		//repulsive force -> constVal * m1 * m2 * lnorm/ lmag*lmag
//			double m1 = _p1.mass, m2 = _p2.mass;
//			myVector lnorm = myVector._normalize(vecL);//unitlength vector of l
//			double fp = -1 * force.constVal * m1 * m2 / (vecL.sqMagn);
//			result[0] = myVector._mult(lnorm,-fp);
//			result[1] = myVector._mult(lnorm,fp);
//		}//only add force if magnitude of distance vector is not 0
//
//		break; }
//	case DAMPSPRING:{//damped spring - not sure if going to use, but what the hey - dependent on old length (need ldot vector)
//		vecL = new myVector(_p2.position.peekFirst(),_p1.position.peekFirst());
//		if (vecL.magn > pa.epsVal) {		//spring with damping force
//			myVector lnorm = myVector._normalize(vecL);//unitlength vector of l
//			myVector lprime = myVector._sub(vecL, myVector._sub(_p1.oldPos.peekFirst(), _p2.oldPos.peekFirst()));		//lprime - time derivative, subtract old length vector from new length vector ?
//			double KsTerm = force.constVal * (vecL.magn - d);
//			double KdTerm = force.constVal2 * (lprime._dot(vecL));
//			double fp = -1 * (KsTerm + KdTerm);
//			result[0] = myVector._mult(lnorm,fp);
//			result[1] = myVector._mult(lnorm, -fp);
//		}//only add force if magnitude of distance vector is not 0
//		break; }
//	case DSPR_THETABAR:{//damped spring to represent ankle
//		vecL = new myVector(_p1.initVel,_p1.velocity.peekFirst());
//		v_l = new myVector(_p2.initPos,_p1.position.peekFirst());
//		if (vecL.magn > pa.epsVal) {		//spring with damping force
//			myVector lnorm = myVector._normalize(vecL);//unitlength vector of l
//			//Eigen::Vector3d lprime = l - (_p1->position[1] - _p2->position[1]);		//lprime - time derivative, subtract old length vector from new length vector ?
//			double KpTerm = force.constVal * (vecL.magn - d);
//			double KdTerm = force.constVal2 * (v_l._dot(vecL));
//			double fp = -1 * (KpTerm + KdTerm);
//			result[0] = myVector._mult(lnorm,fp);
//			result[1] = myVector._mult(lnorm, -fp);
//		}//only add force if magnitude of distance vector is not 0
//		break; }
//	default: {	break; }
//
//}//switch
//return result;

