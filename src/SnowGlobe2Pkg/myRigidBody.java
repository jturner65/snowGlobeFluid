package SnowGlobe2Pkg;

import SnowGlobe2Pkg.SnowGlobeWin.SolverType;

public class myRigidBody{
	protected static mySnowGlobeWin win;
	public final int ID;
	public static int IDgen = 0;

	public myVector initPos, initVel;
	
	public myVector[] aPosition, aLinMomentum, aAngMomentum, aForceAcc, aTorqueAcc, 	
		aOldPos, aOldLinMmnt, aOldAngMmnt, aOldForceAcc, aOldTorqueAcc;
	
    double mMass;
    myQuaternion quatOrient;
    
//	need matrix class
//	Eigen::Matrix3d mOrientation;    
//	Eigen::Matrix3d mMomentofInertia;
//	Eigen::Matrix3d mIBodyInv;
    
	public int curIDX;
	public static final int szAcc = 4;		//size of accumulator arrays
	
	public double mass = 1,origMass;
	public SnowGlobeWin.SolverType solveType;		//{GROUND, EXP_E, MIDPOINT, RK3, RK4, IMP_E, etc }
	public mySolver solver;

	public myRigidBody(mySnowGlobeWin _win, myVector[] _initSt, double _mass, SolverType _styp) {
		win = _win;
		ID = IDgen++;
		mass = _mass;
		init(_initSt, _styp);	
	}
	//init state idxs : 
	//0 : pos, 1 : lin mmnt, 2 : ang mmnt, 3 : force on com, 4 : torque
	private void init(myVector[] _initSt, SnowGlobeWin.SolverType _solv) {
		curIDX = 0;									//cycling ptr to idx in arrays of current sim values
		
		aPosition = new myVector[szAcc];
		aLinMomentum = new myVector[szAcc];
		aAngMomentum = new myVector[szAcc];
		aForceAcc = new myVector[szAcc];
		aTorqueAcc = new myVector[szAcc];
		aOldPos = new myVector[szAcc];
		aOldLinMmnt = new myVector[szAcc];
		aOldAngMmnt = new myVector[szAcc];
		aOldForceAcc = new myVector[szAcc];
		aOldTorqueAcc = new myVector[szAcc];
		
		for(int i=0;i<szAcc;++i){
			aPosition[i] = new myVector();
			aLinMomentum[i] = new myVector();
			aAngMomentum[i] = new myVector();
			aForceAcc[i] = new myVector();
			aTorqueAcc[i] = new myVector();
			aOldPos[i] = new myVector();
			aOldLinMmnt[i] = new myVector();
			aOldAngMmnt[i] = new myVector();
			aOldForceAcc[i] = new myVector();
			aOldTorqueAcc[i] = new myVector();
		}
//		aPosition[0].set(_pos);
//		aLinMomentum[0].set(_velocity);
//		aAngMomentum[0].set(_velocity);
//		aForceAcc[0].set(_forceAcc);
//		aTorqueAcc.set(_forceAcc);
//		aOldPos[0].set(_pos);
//		aOldLinMmnt[0].set(_velocity);
//		aOldAngMmnt[0].set(_velocity);
//		aOldForceAcc[0].set(_forceAcc);
//		aOldTorqueAcc[0].set(_forceAcc);
		
		setOrigMass(mass);
//		initPos = new myVector(_pos);
//		initVel = new myVector(_velocity);
		solveType = _solv;
		solver = new mySolver( _solv);
	}

	private void setOrigMass(double _m) {
		mass = _m;
		origMass = _m;
	}
	public void applyForce(myVector _force) {
		aForceAcc[curIDX]._add(_force);
	}//applyforce
	
	public void intAndAdvance(double deltaT){	
		//TODO rebuild for RB
		
//		myVector[] tSt = new myVector[]{ aPosition[curIDX], aVelocity[curIDX], aOldPos[curIDX], aOldVel[curIDX]};	
//		myVector[] tStDot = new myVector[]{ tSt[1],myVector._div(aForceAcc[curIDX],mass),tSt[3], myVector._div(aOldForceAcc[curIDX],mass)};
//		myVector[] tNSt = solver.Integrate(deltaT, tSt, tStDot);
//		
//		int oldTopIDX = curIDX;
//		curIDX = (curIDX + 1) % szAcc; 
//		aOldPos[curIDX] = aPosition[oldTopIDX];
//		aPosition[curIDX].set(tNSt[0]);
//		
//		aOldVel[curIDX] = aVelocity[oldTopIDX];
//		aVelocity[curIDX].set(tNSt[1]);
//		
//		aOldForceAcc[curIDX] = aForceAcc[oldTopIDX];
//		aForceAcc[curIDX].set(0,0,0);			//clear out new head of force acc	
		
	//	advance(tmpNextStateVec[0], tmpNextStateVec[1], new myVector());
	}
	@Override
	public String toString(){
		String res = "Rigid Body ID : " + ID + "\tMass:"+mass+"\n";
//		res +="\tPosition:"+aPosition[curIDX].toStrBrf()+"\n";
//		res +="\tVelocity:"+aVelocity[curIDX].toStrBrf()+"\n";
//		res +="\tCurrentForces:"+aForceAcc[curIDX].toStrBrf()+"\n";
		return res;		
	}
	
	
	
	
}//myRigidBody

class myBox extends myRigidBody{
	protected static SnowGlobeWin pa;
	public float[] dims;

	public myBox(SnowGlobeWin _pa, mySnowGlobeWin _win, myVector[] _initSt, double _mass, float[] _dims, SolverType _styp) {
		super(_win, _initSt, _mass, _styp);
		pa = _pa;
		dims = _dims;
//        mMomentofInertia(0, 0) = mMass/12.0 * (dims(1) * dims(1) + dims(2) * dims(2));
//        mMomentofInertia(1, 1) = mMass/12.0 * (dims(0) * dims(0) + dims(2) * dims(2));
//        mMomentofInertia(2, 2) = mMass/12.0 * (dims(1) * dims(1) + dims(0) * dims(0));

	}
	
	
	
	
	
}//class myBox