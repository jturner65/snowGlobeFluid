package SnowGlobe2Pkg;

import java.util.*;

import SnowGlobe2Pkg.SnowGlobeWin.SolverType;

public class myParticle {
	public final int ID;
	public static int IDgen = 0;

	public myVector initPos, initVel;
	
	public myVector[] aPosition, aVelocity, aForceAcc, aOldPos, aOldVel, aOldForceAcc;

	public int curIDX;
	public static int curNext, curPrev;
	public int colVal;						//collision check value
	
	public static final int szAcc = 2;		//size of accumulator arrays
	
	public double mass = 1,origMass;
	public SnowGlobeWin.SolverType solveType;		//{GROUND, EXP_E, MIDPOINT, RK3, RK4, IMP_E, etc }
	public mySolver solver;

	public myParticle(myVector _iPos, myVector _iVel, myVector _iFrc, SnowGlobeWin.SolverType _styp) {
		ID = IDgen++;
		init(_iPos, _iVel, _iFrc, _styp);	
	}
	
	private void init(myVector _pos, myVector _velocity, myVector _forceAcc, SnowGlobeWin.SolverType _solv) {
		curIDX = 0;									//cycling ptr to idx in arrays of current sim values
		curNext = 1; 
		curPrev = 0;
		aPosition = new myVector[szAcc];
		aVelocity = new myVector[szAcc];
		aForceAcc = new myVector[szAcc];
		aOldPos = new myVector[szAcc];
		aOldVel = new myVector[szAcc];
		aOldForceAcc = new myVector[szAcc];
		
		for(int i=0;i<szAcc;++i){
			aPosition[i] = new myVector();
			aVelocity[i] = new myVector();
			aForceAcc[i] = new myVector();
			aOldPos[i] = new myVector();
			aOldVel[i] = new myVector();
			aOldForceAcc[i] = new myVector();
		}
		aPosition[0].set(_pos);
		aVelocity[0].set(_velocity);
		aForceAcc[0].set(_forceAcc);
		aOldPos[0].set(_pos);
		aOldVel[0].set(_velocity);
		aOldForceAcc[0].set(_forceAcc);
		
		setOrigMass(mass);
		initPos = new myVector(_pos);
		initVel = new myVector(_velocity);
		solveType = _solv;
		solver = new mySolver( _solv);
	}

	private void setOrigMass(double _m) {
		mass = _m;
		origMass = _m;
	}
	
//	public static void updateCurPtrs(){
//		curNext = (curIDX + 1) % szAcc; 
//		curPrev = 0;
//		
//	}
	
	public void applyForce(myVector _force) {aForceAcc[curIDX]._add(_force);}//applyforce
	public void integAndAdvance(double deltaT){		
		myVector[] tSt = new myVector[]{ aPosition[curIDX], aVelocity[curIDX], aOldPos[curIDX], aOldVel[curIDX]};	
		myVector[] tStDot = new myVector[]{ tSt[1],myVector._div(aForceAcc[curIDX],mass),tSt[3], myVector._div(aOldForceAcc[curIDX],mass)};
		myVector[] tNSt = solver.Integrate(deltaT, tSt, tStDot);
		
		int oldTopIDX = curIDX;
		curIDX = (curIDX + 1) % szAcc; 
		aOldPos[curIDX] = aPosition[oldTopIDX];
		aPosition[curIDX].set(tNSt[0]);
		
		aOldVel[curIDX] = aVelocity[oldTopIDX];
		aVelocity[curIDX].set(tNSt[1]);
		
		aOldForceAcc[curIDX] = aForceAcc[oldTopIDX];
		aForceAcc[curIDX].set(0,0,0);			//clear out new head of force acc	
	}
	
	@Override
	public String toString(){
		String res = "ID : " + ID + "\tMass:"+mass+"\n";
		res +="\tPosition:"+aPosition[curIDX].toStrBrf()+"\n";
		res +="\tVelocity:"+aVelocity[curIDX].toStrBrf()+"\n";
		res +="\tCurrentForces:"+aForceAcc[curIDX].toStrBrf()+"\n";
		
		return res;		
	}
}//myParticle

class mySnowFlake extends myParticle{
	protected static SnowGlobeWin pa;
	protected static mySnowGlobeWin win;
	public int[] color, origColor;	

	public mySnowFlake(SnowGlobeWin _pa, mySnowGlobeWin _win, myVector _iPos, myVector _iVel, myVector _iFrc,
			SolverType _styp) {
		super(_iPos, _iVel, _iFrc, _styp);
		pa = _pa;
		win = _win;
		color = pa.getClr(pa.gui_White, 255);
		origColor = pa.getClr(pa.gui_White, 255);
	}
	
	public void drawMe(){
		pa.pushMatrix();pa.pushStyle();
		pa.translate(aPosition[curIDX]);
		pa.setColorValFill(this.colVal+1, 255);
		pa.sphere(SnowGlobeWin.snowFlakeRad);
		pa.popStyle();pa.popMatrix();
	}
	
	

}//class mySnowFlake