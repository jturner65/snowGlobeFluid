package SnowGlobe2Pkg;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

import SnowGlobe2Pkg.SnowGlobeWin.CollisionType;
import processing.core.*;

public class mySnowGlobeWin extends myDispWindow {
	//ui vars
	public final static int
		gIDX_TimeStep 		= 0,
		gIDX_NumSnowmen 	= 1,
		gIDX_NumStacks		= 2,
		gIDX_NumSnowFlakes  = 3;
	public final int numGUIObjs = 4;												//# of gui objects for ui
	//global variables
	
	//actor values and structures
	public snoMan[] players;                               //our snowmen
	public snoBall[] snoballs;                             //their weapons - 1 per snoman in play at any 1 time
	public stackPresents[] presents;                       //stacks of presents - sets of boxes sitting at random places
	//public int camFollowIdx = -1;                          //index of currently followed snoman - if -1 then centered on middle of screen
	
	//simulation vars - from UI
	private double deltaT = .01;
	public int snoMenCount = 20;                //number of snowmen - min of 3
	public int numStacks = 23;                  //number of present stacks
	private static final int maxNumStacks = 40;
	public int numSnowFlakes = 6000;
	
	public double shakeForceMult = .1; 
	
	//amount of density to add per click/drag 
	public static double addDenseAmt = 10.0;
	
	//snowflakes
	public mySnowFlake[] snowFlakes;
	public myForce[] dfltForces;
	public myCollider[] colliders;
	//interaction value
	public myVector shakeFrcVec;
	public ArrayList<myPoint> shakeDenseVecAra;
	
	//display scaling value
	public final float scVal = 300.0f/SnowGlobeWin.groundRadius;
	
	public snoGlobe globe;									//snow globe construct
	//constants for snowman control and animation
	public final float clickMax = 62831.9f;

	public int snomanFocus=0;                              //which snoman has the spotlight
	public int spotLightChange;                            //counter to change spotlight
	
	public TreeMap<Integer, Float[]> preCalcCos, preCalcSin;
	
	public static final double clickIncr = 0.05;
	
	//to handle real-time update of locations of spheres
	public myVector curMseLookVec;  //pa.c.getMse2DtoMse3DinWorld()
	public myPoint curMseLoc3D;		//pa.c.getMseLoc(pa.sceneCtrVals[pa.sceneIDX])

	//private child-class flags
	public static final int 
			showBndsIDX 	= 0,			//Show fluid bounds values in fluid globe
			showDensIDX 	= 1,			//Show density values in fluid globe
			showVelIDX 		= 2,			//show velocity vectors in fluid globe
			showGndColIDX 	= 3,			//display ground collider TODO
			showSnowFlksIDX = 4,
			showPrsntsIDX 	= 5,
			showSnoMen 		= 6,
			
			isShakeVelIDX 	= 7,			
			isShakeDenseIDX = 8;
	public static final int numPrivFlags = 9;
	
	
	public mySnowGlobeWin(SnowGlobeWin _p, String _n, int _flagIdx, int[] fc, int[] sc, float[] rd, float[] rdClosed,String _winTxt, boolean _canDrawTraj) {
		super(_p, _n, _flagIdx, fc, sc, rd, rdClosed, _winTxt, _canDrawTraj);
		float stY = rectDim[1]+rectDim[3]-4*yOff,stYFlags = stY + 2*yOff;
		trajFillClrCnst = SnowGlobeWin.gui_DarkCyan;		//override this in the ctor of the instancing window class
		trajStrkClrCnst = SnowGlobeWin.gui_Cyan;
		super.initThisWin(_canDrawTraj, true);
	}
	
	@Override
	//initialize all private-flag based UI buttons here - called by base class
	public void initAllPrivBtns(){
		truePrivFlagNames = new String[]{								//needs to be in order of flags
				"Showing Fluid Bounds","Showing Fluid Densities","Showing Fluid Velocities","Showing Ground Collider",
				"Showing Snowflakes","Showing Presents","Showing Snowmen" 
		};
		falsePrivFlagNames = new String[]{			//needs to be in order of flags
				"Hiding Fluid Bounds","Hiding Fluid Densities","Hiding Fluid Velocities","Hiding Ground Collider",
				"Hiding Snowflakes","Hiding Presents","Hiding Snowmen"
		};
		privModFlgIdxs = new int[]{showBndsIDX,showDensIDX,showVelIDX,showGndColIDX,showSnowFlksIDX,showPrsntsIDX,showSnoMen};
		numClickBools = privModFlgIdxs.length;	
		//maybe have call for 		initPrivBtnRects(0);		 here
		initPrivBtnRects(0,numClickBools);
	}//initAllPrivBtns
	
	@Override
	protected void initMe() {	
		initUIBox();				//set up ui click region to be in sidebar menu below menu's entries		
		dispFlags[trajDecays] = true;								//this window responds to travelling reticle/playing
		curTrajAraIDX = 0;
		initPrivFlags(numPrivFlags);
		initAllSnowGlobe();
		//this.snomanFocus = (int)(Math.random()* snoMenCount);    //randomly pick a snoman to get a spotlight 
		this.spotLightChange = 0; 
		
		setPrivFlags(showSnowFlksIDX, true);		
		setPrivFlags(showPrsntsIDX, true);		
		setPrivFlags(showSnoMen, true);	
	}
	
	//initialize all snowglobe-specific variables
	public void initAllSnowGlobe(){
		//initialize the array of snomen	
		shakeFrcVec = new myVector();		//direction being shaken
		shakeDenseVecAra = new ArrayList<myPoint>();		//location for density to be increased
		globe = new snoGlobe(pa,this, pa.globeCtr);
		initSnowmen();		
		buildDefForces("snowglobe",-4.0);
		buildDefColliders();
		initSnowflakes();
		initPresents();
	}//init snomen
	
	private void initPresents(){
		presents = new stackPresents[maxNumStacks];		
		//initialize stacks of presents
		for (int i = 0; i < maxNumStacks; ++i){  
			TreeMap<Float, stackPresents> tmpStack = new TreeMap<Float, stackPresents>(),
					ptnStack = new TreeMap<Float, stackPresents>();
			boolean tooClose = true;
			int iter = 0;
			stackPresents tmp = null;
			while((iter < 100) && (tooClose)){
				tmp = new stackPresents(pa,i);
				tooClose = false;
				tmpStack = new TreeMap<Float, stackPresents>();
				for(int j=0;j<i;++j){tmpStack.put(myVectorf._dist(tmp.location, presents[j].location), presents[j]);}
				if((tmpStack.size() > 0) && (tmpStack.firstKey() < (tmp.getStackWidth() + tmpStack.firstEntry().getValue().getStackWidth()))){
					ptnStack.put(tmpStack.firstKey(),tmp);//closest distance to neighbor - use this to get stack furthest from closest neighbor if thresholding is too strict
					tooClose = true;
					//pa.pr("Too close : " + tmpStack.firstKey()+ " between : " + tmp.ID+" and " + tmpStack.firstEntry().getValue().ID);
				}
				iter++;
			}//while			
			presents[i] = (tmp != null ? tmp : ptnStack.lastEntry().getValue());//if null get package that was furthest in initial check
		}
	}
	
	//call after numstacks value changed
	private void reinitPresents(){
		for (int i = 0; i < numStacks; ++i){  
			TreeMap<Float, stackPresents> tmpStack = new TreeMap<Float, stackPresents>(),
					ptnStack = new TreeMap<Float, stackPresents>();
			boolean tooClose = true;
			int iter = 0;
			while((iter < 100) && (tooClose)){
				presents[i].resetLocation();
				tooClose = false;
				tmpStack = new TreeMap<Float, stackPresents>();
				for(int j=0;j<i;++j){tmpStack.put(myVectorf._dist(presents[i].location, presents[j].location), presents[j]);}
				if((tmpStack.size() > 0) && (tmpStack.firstKey() < (presents[i].getStackWidth() + tmpStack.firstEntry().getValue().getStackWidth()))){
					ptnStack.put(tmpStack.firstKey(),presents[i]);//closest distance to neighbor - use this to get stack furthest from closest neighbor if thresholding is too strict
					tooClose = true;
					//pa.pr("Too close : " + tmpStack.firstKey()+ " between : " + tmp.ID+" and " + tmpStack.firstEntry().getValue().ID);
				}
				iter++;
			}//while			
		}
	}
	
	private void initSnowflakes(){
		snowFlakes = new mySnowFlake[numSnowFlakes];
		for(int i=0;i<numSnowFlakes;++i){
			//(SnowGlobeWin _pa, mySnowGlobeWin _win, myVector _iPos, myVector _iVel, myVector _iFrc, SnowGlobeWin.SolverType _styp)
			//snowFlakes[i] = new myParticle(pa,this, getRandPos(snoGlobe.snowGlobRad), new myVector() , new myVector(), SnowGlobeWin.SolverType.RK4);
			snowFlakes[i] = new mySnowFlake(pa,this, getRandPos(snoGlobe.snowGlobRad), new myVector() , new myVector(), SnowGlobeWin.SolverType.RK4);
		}
	}
	
	private void initSnowmen(){
		players = new snoMan[snoMenCount];
		snoballs = new snoBall[snoMenCount];
		//set up a center vector for where the snoman will be centered
		myVectorf smCenter = new myVectorf(0,0,0); 
		float smHeight;
		//generate starting values for the snomen
		//and make each one 
		for (int i = 0; i < snoMenCount; ++i){
			smHeight = (float) (pa.avgSMHeight + ((Math.random()*(pa.avgSMHeight/2.0)) - (pa.avgSMHeight/4.0)));
			smCenter.set(.7*((Math.random()*( SnowGlobeWin.groundRadius*2)) -  SnowGlobeWin.groundRadius),.7*((Math.random()*( SnowGlobeWin.groundRadius*2)) -  SnowGlobeWin.groundRadius),0);                   //set center of this snowman's butt
			float thet = pa.random(PConstants.TWO_PI);
			players[i] = new snoMan(pa,this, smCenter, smHeight, thet, i);
			snoballs[i] = new snoBall(pa,this,i);
		}//for int i - creates array of snowmen		
	}	
	
	//random position within snowglobe - z is vertical
	private static final double third = 1.0/3.0;
	private myVector getRandPos(double rad){
		myVector pos = new myVector();
		do{
			double u = ThreadLocalRandom.current().nextDouble(0,1), r = rad * Math.pow(u, third),
					cosTheta = ThreadLocalRandom.current().nextDouble(-1,1), sinTheta =  Math.sin(Math.acos(cosTheta)),
					phi = ThreadLocalRandom.current().nextDouble(0,PConstants.TWO_PI);
			pos.set(sinTheta * Math.cos(phi), sinTheta * Math.sin(phi),cosTheta);
			pos._mult(r);
			pos._add(globe.center);
		} while (pos.z < 0);
		return pos;
	}
	
	@Override
	//set flag values and execute special functionality for this sequencer
	public void setPrivFlags(int idx, boolean val){
		int flIDX = idx/32, mask = 1<<(idx%32);
		privFlags[flIDX] = (val ?  privFlags[flIDX] | mask : privFlags[flIDX] & ~mask);
		switch(idx){
			case showBndsIDX : {
				pa.pr("showBndsIDX set : " + val);
				break;}
			case showDensIDX : {
				pa.pr("showDensIDX set : " + val);
				break;}				//a sphere has been selected and either fixed or released
			case showVelIDX : {
				pa.pr("showVelIDX set : " + val);
				break;}
			case showGndColIDX : {
				pa.pr("showGndColIDX set : " + val);
				break;}
			case showSnowFlksIDX : {
				pa.pr("showSnowFlksIDX set : " + val);
				break;}
			case showPrsntsIDX : {
				pa.pr("showPrsntsIDX set : " + val);
				break;}
			case showSnoMen : {
				pa.pr("showSnoMen set : " + val);
				break;}
			case isShakeVelIDX : {
				//pa.pr("isShakeVelIDX set : " + val);
				break;}
			case isShakeDenseIDX : {
				//pa.pr("isShakeDenseIDX set : " + val);
				break;}
		}		
	}//setPrivFlags		
	
	//initialize structure to hold modifiable menu regions
	@Override
	protected void setupGUIObjsAras(){	
		//pa.outStr2Scr("setupGUIObjsAras in :"+ name);
		guiMinMaxModVals = new double [][]{//min max mod values for each modifiable UI comp
			{.001, .1, .001},								//timestep
			{1, 100, 1},									//# of snowmen
			{1, maxNumStacks, 1},										//# of presents
			{100,10000,100}									//# of snowflakes
		};						
		guiStVals = new double[]{								//starting value
				.01,
				20,
				23,
				6000
		};
		
		guiObjNames = new String[]{							//name/label of component
				"Time Step",
				"# of Snowmen",
				"# of Presents Stacks",
				"# of Snowflakes"
		};
		//idx 0 is treat as int, idx 1 is obj has list vals, idx 2 is object gets sent to windows
		guiBoolVals = new boolean [][]{
			{false, false, true},
			{true, false, true},
			{true, false, true},
			{true, false, true}
			};						//per-object  list of boolean flags
		
		//since horizontal row of UI comps, uiClkCoords[2] will be set in buildGUIObjs		
		guiObjs = new myGUIObj[numGUIObjs];			//list of modifiable gui objects
		if(numGUIObjs > 0){
			buildGUIObjs(guiObjNames,guiStVals,guiMinMaxModVals,guiBoolVals,new double[]{xOff,yOff});			//builds a horizontal list of UI comps
		}
	}
	
	@Override
	protected void setUIWinVals(int UIidx) {
		double val = guiObjs[UIidx].getVal();
		int ival = (int)val;
		switch(UIidx){
		case gIDX_TimeStep			: {
			if(val != deltaT){deltaT = val;}
			break;}
		case gIDX_NumSnowmen 	: {
			if(ival != snoMenCount){
				snoMenCount = ival;
				pa.setFlags(pa.runSim, false);
				initSnowmen();
			}
			break;}
		case gIDX_NumStacks 	: {		
			if(ival != numStacks){
				numStacks = ival;
				pa.setFlags(pa.runSim, false);
				reinitPresents();
			}
			break;}
		case gIDX_NumSnowFlakes : {
			if(ival != numSnowFlakes){
				numSnowFlakes = ival;
				pa.setFlags(pa.runSim, false);
				initSnowflakes();
			}
			break;}
		default : {break;}
		}
	}//setUIWinVals
	
	//if any ui values have a string behind them for display
	@Override
	protected String getUIListValStr(int UIidx, int validx) {
		switch(UIidx){
			default : {break;}
		}
		return "";
	}


	@Override
	public void initDrwnTrajIndiv(){
	}

	//overrides function in base class mseClkDisp
	@Override
	public void drawTraj3D(float animTimeMod,myPoint trans){
		//pa.outStr2Scr("mySphereUI.drawTraj3D() : I am overridden in 3d instancing class", true);
	}//drawTraj3D

	@Override
	protected void drawMe(float animTimeMod) {
		pa.background(0,15,55,255);
		curMseLookVec = pa.c.getMse2DtoMse3DinWorld(pa.sceneCtrVals[pa.sceneIDX]);			
		curMseLoc3D = pa.c.getMseLoc(pa.sceneCtrVals[pa.sceneIDX]);
		pa.pushMatrix();pa.pushStyle();
		pa.scale(scVal);
		pa.noLights();											
		drawLights(animTimeMod*.5f , snomanFocus);                   //make lighting
		pa.pushMatrix();pa.pushStyle();
		pa.noStroke();
		if(getPrivFlags(showSnowFlksIDX)){drawSnowFlakes(); }		
		if(pa.flags[pa.debugMode]){		globe.drawSnowGlobeDebug();} 
		else {							
			if(getPrivFlags(showSnoMen)){drawSnoMen();}
			if(getPrivFlags(showPrsntsIDX)){drawPresents();}//for the stacks that are more than the count of the snomen	
			globe.drawSnowGlobe();
		}
		pa.popStyle();pa.popMatrix();
		pa.noLights();
		//pa.turnOnLights();		
		if(pa.flags[pa.runSim]){executeSnowGlobe();}				//need to build better mechanism to be draw-driven instead of 1-shot for play
		//pa.pr("win draw");
		pa.popStyle();pa.popMatrix();
	}
	
	public void executeSnowGlobe(){//runtime simulation stuff 
		//apply forces
		applyForcesToSystem();
		
		if (getPrivFlags(isShakeVelIDX)) {
			setPrivFlags(isShakeVelIDX, false);
			//addShakeForceToFluid(fMult, msdrgVal0, msdrgVal1); 		getPrivFlags(this.isShakenIDX] = false;			//setting false keeps force from compounding
			globe.myFluidBoxAddForce(shakeFrcVec);
			shakeFrcVec.set(0,0,0);
		}
		if(getPrivFlags(isShakeDenseIDX)){
			setPrivFlags(isShakeDenseIDX, false);
			for (myPoint vec : shakeDenseVecAra){				
				globe.myFluidBoxAddDens( addDenseAmt, vec);
			}
			shakeDenseVecAra.clear();
		}

		globe.timeStep();
		
		calcFluidForceForAllParticles();
		
		handlePartCldrCollision();
		for (snoMan s : players){	s.simulation();}		
		invokeSolverDerivEval();		
	}//executeSnowGlobe

	//check all collisions - particle to particle and square/spherical boundary(ec)
	public void handlePartCldrCollision(){//check for each collider for all particles - check returns 0 - no collision possible, 1- breach, 2- collision due to velocity, 3- collision due to contact force
       // myVector resFrc = new myVector(), tmpFrc;
	
		for (int partIDX = 0; partIDX < snowFlakes.length; ++partIDX) {
			for (int i = 0; i < colliders.length; ++i) {
				snowFlakes[partIDX].colVal = colliders[i].checkCollision(deltaT, snowFlakes[partIDX]);
				if(snowFlakes[partIDX].colVal > 0){	
					colliders[i].handleCollision(snowFlakes[partIDX], snowFlakes[partIDX].colVal);	
					continue;
				}
			}//for each particle
		}//for each collider		
 	}//check for collisions
	
	private final double log11 = Math.log(1.1), logOffset = 6;
	//calculate the effects of the fluid sim on the particles in the globe - to take the place of rigged forces in shaking routine
	public void calcFluidForceForAllParticles() {
		myVector fluidFrcVal = new myVector(), applyShakeVal, tmpFrc1;
		double //shakeValMag,
			logShakeValp1;
		int curPartIDX = snowFlakes[0].curIDX;
		for (int pidx = 0; pidx < snowFlakes.length; ++pidx) {
			fluidFrcVal = globe.getVelAtCell(snowFlakes[pidx].aPosition[curPartIDX]);
			if (fluidFrcVal.magn > .00001) {
				logShakeValp1 = 1 + Math.log(fluidFrcVal.magn) / log11;
				if(logShakeValp1 > logOffset){
					applyShakeVal = myVector._mult(myVector._normalize(fluidFrcVal), logOffset + Math.log(logShakeValp1)) ;					
				} else {
					applyShakeVal = myVector._mult(myVector._normalize(fluidFrcVal), Math.abs (logShakeValp1));					
				}
				snowFlakes[pidx].applyForce(applyShakeVal);// * (colliders[0].getDistFromCenter(p[pidx]->position[0], 1)/(1.0 * (snowGlobRad))));
			}
		}//for all particles
	}//calcFluidForceForAllParticles

	//default forces that always affect particles
	public void buildDefForces(String name, double kdrag) {
		myForce gravity = new mySclrForce(pa, this, name + "Force_Grav", SnowGlobeWin.gravVec);
		if (kdrag != 0) {
			dfltForces = new myForce[]{gravity, new myVecForce(pa, this, name + "Fluid_Drag", kdrag)};
		} else {
			dfltForces = new myForce[]{gravity};			
		}
	}//default forces
	
	//build default colliders for snow globe and ground in snowglobe - build after globe
	public void buildDefColliders(){
		colliders = new myCollider[2];//ground and sphere - TODO handle collisions with packages
		float vertCoord = snoGlobe.snowGlobRad * 10.0f;
//		systems[curSystemIDX]->buildGndCollider(.00001, 0, 10, gndLocY, 10, std::string("GroundPart1"), Eigen::Vector3d(0, gndLocY, -12), Eigen::Vector3d(0, gndLocY, -12));
//		systems[curSystemIDX]->buildGlobeCollider(.00001, 0, snowGlobRad, distFromSnowGlobe);
		
		//buildGndCollider(double kRest, double muFrict, double x, double y, double z, String name,myVector drLoc, myVector gndLoc) {
		myVector[] vertLocs = new myVector[]{
				new myVector(vertCoord, vertCoord, 0),
				new myVector(-vertCoord, vertCoord, 0),
				new myVector(-vertCoord, -vertCoord, 0),
				new myVector(vertCoord, -vertCoord, 0)};
		//ground
		colliders[0] = new planeCollider(pa, this, name, 
				new myVectorf(globe.center), 													//draw location - doesn't  matter - should draw collider through passed points
				vertLocs); 												//verticies on plane - can be any values so long as height is correct
		colliders[0].Krest = .00001;
		colliders[0].muFrict = 0;        //friction force
			
		//buildGlobeCollider(double krest, double muFrict, double rad, double distFromGlb) 
		colliders[1] = new sphereCollider(pa, this,"SnowGlobePart2", 
				new myVectorf(globe.center),								//drawLoc
				new myVector(globe.center),				//center
				new myVector(snoGlobe.snowGlobRad, snoGlobe.snowGlobRad, snoGlobe.snowGlobRad), 								//radius
				true);	
		colliders[1].Krest = .00001;
		colliders[1].muFrict = 0;
	}
	
//	public void addShakeForceToFluid(double fMult, myVector msdrgVal0, myVector msdrgVal1) {
//		shakeVal = myVector._mult(myVector._sub(msdrgVal1, msdrgVal0), fMult);
//		globe.myFluidBoxAddForce(shakeVal);
//	}//addForceToFluid		
	
	//this will apply all system-wide forces to all particles - individual forces set by alternate methods
	public void applyForcesToSystem() {
		myVector[] result;
		for (int fidx = 0; fidx < dfltForces.length; ++fidx) {
			//if(null == dfltForces[fidx]){continue;}
			for (int pidx = 0; pidx < snowFlakes.length; ++pidx) {
				//second particle and d ignored for single particle forces
				result = dfltForces[fidx].calcForceOnParticle(snowFlakes[pidx], snowFlakes[pidx], 0);//do not use d for single particle forces
				//idx 1 is force on other particle - if we have multi particle forces we might be able to minimize overhead by calculating only 1 time for 2 parts
				snowFlakes[pidx].applyForce(result[0]);
			}//for each force
		}//for all particles
	}//applyForcesToSystem
	
	//call after fluid and collision forces
	public void invokeSolverDerivEval() {for (int idx = 0; idx < snowFlakes.length; ++idx) {snowFlakes[idx].intAndAdvance(deltaT);}}
	
	private void stepSpotLight(){
		if(spotLightChange % 50 == 0){                    //change spotlight after a little while
			spotLightChange = 0;
			snomanFocus = (snomanFocus+1)%2;                //go from even to odd by toggling between 1 and 0
		}
		spotLightChange++;
	}//stepSpotLight
	
	public void drawLights(double click, int smIdxArg){
		stepSpotLight();
		pa.lightSpecular(255, 255, 255);     
		//spotLight(v1, v2, v3, x, y, z, nx, ny, nz, angle, concentration)
		myVector lightDir = new myVector(players[smIdxArg].location.x, players[smIdxArg].location.y,players[smIdxArg].location.z-500);
		pa.spotLight(250,250,250,0,0,500,(float)lightDir.x, (float)lightDir.y, (float)lightDir.z, PConstants.PI/5.0f, 1000);
		pa.ambientLight(110,110,110);
		//spinning light - need a disco ball :D
		pa.directionalLight(255, 0, 0, (float)Math.sin(click), (float)Math.cos(click), -1);
		pa.directionalLight(0, 0, 255, (float) Math.sin(click + PConstants.HALF_PI),(float) Math.cos(click + PConstants.HALF_PI), -1);
		pa.directionalLight(0, 255, 0, (float)(-1*Math.sin(click - PConstants.HALF_PI)), (float)(-1*Math.cos(click - PConstants.HALF_PI)), -1);
	}//createLight
		
	private void drawSnowFlakes(){
		pa.pushMatrix();pa.pushStyle();
		pa.sphereDetail(4);
		for(int i=0;i<numSnowFlakes;++i){	snowFlakes[i].drawMe();	}
		pa.popStyle();pa.popMatrix();
	}//drawSnowFlakes
	
	//draw the snowmen, during draw routine
	public void drawSnoMen(){
		pa.pushMatrix();pa.pushStyle();
		//for (int i = 0; i < snoMenCount; ++i){	players[i].drawMe();	}		
		for (snoMan s : players){	s.drawMe();	}		
		pa.popStyle();pa.popMatrix();
	}//drawSnoMen
	
	
	//draw the presents, during draw routine
	public void drawPresents(){
		pa.pushMatrix();pa.pushStyle();
		for (int i = 0; i < numStacks; ++i){ presents[i].drawMe();   	}//for the stacks that are more than the count of the snomen				
		pa.popStyle();pa.popMatrix();
	}//drawSnoMen
	
	/**
	* checks if passed snoman is colliding with any other snomen in the scene - need to modify to build avoidance vector based on dist to all other snowmen
	* @param location location vector to check
	* @return the snoman that this is colliding with
	*/
	public int collidingSM(snoMan mover){
		myVectorf link = new myVectorf();
		float moveR = pa.groundRadius/50.0f;
		for(int i = 0; i < snoMenCount; ++i){
			if (mover.ID == i){continue;}
			link.set(players[i].location,mover.location);
			if (link.magn < (mover.buttRadL + players[i].buttRadL)){
				link._scale((float)Math.random()*moveR);//
				mover.location._add(link);
				mover.heading.set(link);
				mover.heading._normalize();
				//collision happened - displace some direction away from players[i] - vector is from players[i] to mover

				return i;
			}//check if colision  
		}//for i < pa.snoMenCount  
		return -1; 
	}//colliding snomen method
	
	
	@Override
	protected void playMe() {	}//only called 1 time
	//when animation stops
	@Override
	protected void stopMe() {}	
	
	//print out all trajectories for current sphere
	public void dbgFunction1(){

	}
	
	public void dbgFunction2(){		

	}
	
	@Override
	public void clickDebug(int btnNum){
		pa.outStr2Scr("click debug in "+name+" : btn : " + btnNum);
		switch(btnNum){
			case 0 : {//note vals
				dbgFunction1();
				break;
			}
			case 1 : {//measure vals
				dbgFunction2();				
				break;
			}
			case 2 : {
				break;
			}
			case 3 : {
				break;
			}
			default : {break;}
		}		
	}
	@Override
	protected void processTrajIndiv(myDrawnSmplTraj drawnTraj){
		pa.outStr2Scr("Process traj in snow globe");
		//traj processing
	}
	@Override
	protected myPoint getMouseLoc3D(int mouseX, int mouseY) {
		return pa.P(-100000,10000,100000);			//dummy point - remove this crap eventually when handle trajs correctly
	}//getMouseLoc3D
	@Override
	protected boolean hndlMouseMoveIndiv(int mouseX, int mouseY, myPoint mseClckInWorld){
		return false;
	}
	//alt key pressed handles trajectory
	//cntl key pressed handles unfocus of spherey
	@Override
	//mseBtn : 0 is left, 1 is right button
	protected boolean hndlMouseClickIndiv(int mouseX, int mouseY, myPoint mseClckInWorld, int mseBtn) {
		boolean res = checkUIButtons(mouseX, mouseY);
		if(res) {return res;}
		//pa.outStr2Scr("sphere ui click in world : " + mseClckInWorld.toStrBrf());
		if(mseBtn == 1){			
			shakeDenseVecAra.add(new myPoint(mseClckInWorld)); //add to velocity
			if(!getPrivFlags(isShakeDenseIDX)) {setPrivFlags(isShakeDenseIDX, true);}		
			res = true;
		}//add to density at click location
		return res;
	}//hndlMouseClickIndiv

	@Override
	//mseBtn : 0 is left, 1 is right button
	protected boolean hndlMouseDragIndiv(int mouseX, int mouseY, int pmouseX, int pmouseY, myPoint mouseClickIn3D, myVector mseDragInWorld, int mseBtn) {
		boolean res = false;
		//pa.pr("win mouse drag");
		if(mseBtn == 0){	
			myVector tmp = myVector._mult(mseDragInWorld, shakeForceMult); //add to velocity
			if(tmp.magn > pa.epsVal){shakeFrcVec = tmp;}
			if(!getPrivFlags(isShakeVelIDX)){setPrivFlags(isShakeVelIDX, true);}
			res = true;
		} else if(mseBtn == 1) {					 
			shakeDenseVecAra.add(mouseClickIn3D); //add to velocity
			if(!getPrivFlags(isShakeDenseIDX)){setPrivFlags(isShakeDenseIDX, true);}
			res = true;
		}//add to density
		return res;
	}//hndlMouseDragIndiv
	
	@Override
	protected void snapMouseLocs(int oldMouseX, int oldMouseY, int[] newMouseLoc) {}	

	@Override
	protected void hndlMouseRelIndiv() {}
	@Override
	protected void endShiftKeyI() {}
	@Override
	protected void endAltKeyI() {}
	@Override
	protected void endCntlKeyI() {}
	@Override
	protected void addSScrToWinIndiv(int newWinKey){}
	@Override
	protected void addTrajToScrIndiv(int subScrKey, String newTrajKey){}
	@Override
	protected void delSScrToWinIndiv(int idx) {}	
	@Override
	protected void delTrajToScrIndiv(int subScrKey, String newTrajKey) {}
	//resize drawn all trajectories
	@Override
	protected void resizeMe(float scale) {}
	@Override
	protected void closeMe() {}
	@Override
	protected void showMe() {}
}
