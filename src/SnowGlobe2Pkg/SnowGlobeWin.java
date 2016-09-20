package SnowGlobe2Pkg;

import java.awt.event.KeyEvent;
import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;			//used for threading
import java.util.concurrent.Executors;

import processing.core.*;
import processing.opengl.*;
/**
 * Snowglobe with snowmen
 * John Turner
 */


 public class SnowGlobeWin extends PApplet {
	 
	public enum ForceType {
		F_NONE(0), S_SCALAR(1), S_VECTOR(2), ATTR(3), REPL(4), DAMPSPRING(5), DSPR_THETABAR(6);		
		private int value; 
		private static Map<Integer, ForceType> map = new HashMap<Integer, ForceType>(); 
	    static { for (ForceType enumV : ForceType.values()) { map.put(enumV.value, enumV);}}
		private ForceType(int _val){value = _val;} 
		public int getVal(){return value;} 	
		public static ForceType getVal(int idx){return map.get(idx);}
		public static int getNumVals(){return map.size();}						//get # of values in enum	
	};
	public enum ConstraintType {
		C_NONE(0), C_Circular(1);
		private int value; 
		private static Map<Integer, ConstraintType> map = new HashMap<Integer, ConstraintType>(); 
	    static { for (ConstraintType enumV : ConstraintType.values()) { map.put(enumV.value, enumV);}}
		private ConstraintType(int _val){value = _val;} 
		public int getVal(){return value;} 	
		public static ConstraintType getVal(int idx){return map.get(idx);}
		public static int getNumVals(){return map.size();}
	};
	
	public enum CollisionType {
		CL_NONE(0), FLAT(1), PARTICLE(2), SPHERE(3), CYLINDER(4), BOX(5);
		private int value; 
		private static Map<Integer, CollisionType> map = new HashMap<Integer, CollisionType>(); 
	    static { for (CollisionType enumV : CollisionType.values()) { map.put(enumV.value, enumV);}}
		private CollisionType(int _val){value = _val;} 
		public int getVal(){return value;} 	
		public static CollisionType getVal(int idx){return map.get(idx);}
		public static int getNumVals(){return map.size();}						//get # of values in enum	
	};

	public enum SolverType {
		GROUND(0), EXP_E(1), MIDPOINT(2), RK3(3), RK4(4), IMP_E(5), TRAP(6), VERLET(7), RK4_G(8);
		private int value; 
		private static Map<Integer, SolverType> map = new HashMap<Integer, SolverType>(); 
	    static { for (SolverType enumV : SolverType.values()) { map.put(enumV.value, enumV);}}
		private SolverType(int _val){value = _val;} 
		public int getVal(){return value;} 	
		public static SolverType getVal(int idx){return map.get(idx);}
		public static int getNumVals(){return map.size();}						//get # of values in enum			
	};
	
	public static String[] ForceType2str = {"None", "Gravity-type force (scalar particle quantity)", "Air Drag-type force (vector particle quantity)", "Attractor", "Repulsor", "Damped Spring", "Force back to ThetaBar"};
	public static String[] ConstraintType2str = {"None", "Circular/Bar-type constraint"};
	public static String[] CollisionType2str = {"None", "Flat surface", "Particle to particle", "Spherical", "Cylinder", "Bounding Box"};
	public static String[] SolverType2str = {"Ground_Truth", "Explicit_Euler", "Midpoint", "RK3", "RK4", "Implicit_Euler", "Trapezoidal", "Verlet", "Gen_RK4"};

    public static myVector gravVec = new myVector(0, 0, -9.8);		//derp z is up/down in globe world

	 
	//project-specific variables
	public String prjNmLong = "Snow Globe Advanced with JCuda", prjNmShrt = "SnowGlobeJCuda";
	
	//global variables		
	public final int drawnTrajEditWidth = 10; //TODO make ui component			//width in cntl points of the amount of the drawn trajectory deformed by dragging
	public final float
				wScale = frameRate/5.0f,					//velocity drag scaling	
				trajDragScaleAmt = 100.0f;					//amt of displacement when dragging drawn trajectory to edit
			
	public String msClkStr = "";
			
	public int glblStartSimTime, glblLastSimTime;
	
	//constants
	public static final float groundRadius = 30.0f;           //radius of the disk of the ground - all components need to be a multiple of this
	public static myVectorf globeCtr = new myVectorf(0,0, SnowGlobeWin.groundRadius/2.0);

	
	public static final double timeIncr = .05;               //amount time increments per cycle of draws
//	public static final int numStacks = 23;                  //number of present stacks
//	public static final int snoMenCount = 20;                //number of snowmen - min of 3
	public static final float snowFlakeRad = groundRadius * .004f, tenSFRads = 10*snowFlakeRad;
	
	public static final float gSqRad = groundRadius*.95f * groundRadius *.95f,
				avgSMHeight = SnowGlobeWin.groundRadius * .10667f;
		
	
/////////////////////////////		
//CODE STARTS
////////////////////////////
	public void settings(){	size((int)(displayWidth*.95f), (int)(displayHeight*.9f),P3D);	}
	public void setup(){
		initOnce();
		background(bground[0],bground[1],bground[2],bground[3]);		
	}//setup	
	
	//called once at start of program
	public void initOnce(){
		initVisOnce();						//always first
		sceneIDX = 1;//(flags[show3D] ? 1 : 0);
		glblStartSimTime = millis();
		glblLastSimTime =  millis();
		
		numThreadsAvail = Runtime.getRuntime().availableProcessors();
		pr("# threads : "+ numThreadsAvail);
		th_exec = Executors.newCachedThreadPool();
		
		focusTar = new myVector(sceneFcsVals[sceneIDX]);

		initDispWins();
		setFlags(showUIMenu, true);					//show input UI menu	
		setFlags(showSnowGlbWin, true);				//show snowglobe sim
		setCamView(); 
		initProgram();
	}//initOnce
	
	//called multiple times, whenever re-initing
	public void initProgram(){
		initVisProg();				//always first
		drawCount = 0;
	}//initProgram

	public void draw(){	
		animCntr = (animCntr + (baseAnimSpd )*animModMult) % maxAnimCntr;						//set animcntr - used only to animate visuals		
		//cyclModCmp = (drawCount % ((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiObjs[((mySideBarMenu)dispWinFrames[dispMenuIDX]).gIDX_cycModDraw].valAsInt() == 0);
		pushMatrix();pushStyle();
		drawSetup();																//initialize camera, lights and scene orientation and set up eye movement
//		//if ((!cyclModCmp) || (flags[runSim])) {drawCount++;}						//needed to stop draw update so that pausing sim retains animation positions			
//		if ((flags[runSim])) {drawCount++;}											//needed to stop draw update so that pausing sim retains animation positions			
//		glblStartSimTime = millis();
//		//float modAmtSec = (glblStartSimTime - glblLastSimTime)/1000.0f;
//		glblLastSimTime = millis();
////		if(flags[runSim] ){
////			//run simulation
////		}		//play in current window
		if((curFocusWin == -1) || (dispWinIs3D[curFocusWin])){		
			background(bground[0],bground[1],bground[2],bground[3]);				//if refreshing screen, this clears screen, sets background
			translate(focusTar.x,focusTar.y,focusTar.z);								//focus location center of screen		
			draw3D();
		} 		//if not 3d scene, draw drawUI
		buildCanvas();		
		popStyle();popMatrix(); 
		drawUIAnd2D();																	//draw UI overlay on top of rendered results			
		if (flags[saveAnim]) {	savePic();}
		consoleStrings.clear();
		surface.setTitle(prjNmLong + " : " + (int)(frameRate) + " fps|cyc \t sceneIDX : " + sceneIDX);
	}//draw
	
	public void draw3D(){
		pushMatrix();pushStyle();
		translateSceneCtr();				//move to center of 3d volume to start drawing	
//		for(int i =1; i<numDispWins; ++i){
		if((isShowingWindow(1)) && (dispWinFrames[1].dispFlags[myDispWindow.is3DWin])){dispWinFrames[1].draw(myPoint._add(sceneCtrVals[sceneIDX],focusTar));}
//		}
		popStyle();popMatrix();
		drawAxes(100,3, new myPoint(-c.viewDimW/2.0f+40,0.0f,0.0f), 200, false); 		//for visualisation purposes and to show movement and location in otherwise empty scene
		if(canShow3DBox[this.curFocusWin]) {drawBoxBnds();}
	}//draw3D	

	public void buildCanvas(){
		c.buildCanvas();
		c.drawMseEdge();
	}
	
	//if should show problem # i
	public boolean isShowingWindow(int i){return flags[(i+this.showUIMenu)];}//showUIMenu is first flag of window showing flags
	public void drawUIAnd2D(){					
		for(int i =1; i<numDispWins; ++i){if ( !(dispWinFrames[i].dispFlags[myDispWindow.is3DWin])){dispWinFrames[i].draw(sceneCtrVals[sceneIDX]);}}
		//dispWinFrames[0].draw(sceneCtrVals[sceneIDX]);
		for(int i =1; i<numDispWins; ++i){dispWinFrames[i].drawHeader();}
		//menu
		dispWinFrames[0].draw(sceneCtrVals[sceneIDX]);
		dispWinFrames[0].drawHeader();
		drawOnScreenData();				//debug and on-screen data
		//noLights();
		displayHeader();				//my pic and name
		//lights();
	}//drawUI	
	public void translateSceneCtr(){translate(sceneCtrVals[sceneIDX].x,sceneCtrVals[sceneIDX].y,sceneCtrVals[sceneIDX].z);}
	
//		//perform full translations for 3d render - scene center and focus tar
//		public void translateFull3D(){
//			translate(focusTar.x,focusTar.y,focusTar.z);
//			translateSceneCtr();
//		}
	
	public void setFocus(){
		focusTar.set(sceneFcsVals[(sceneIDX+sceneFcsVals.length)%sceneFcsVals.length]);
		switch (sceneIDX){//special handling for each view
		case 0 : {initProgram();break;} //refocus camera on center
		case 1 : {initProgram();break;}  
		}
	}
	
	public void setCamView(){//also sets idx in scene focus and center arrays
		sceneIDX = (curFocusWin == -1 || this.dispWinIs3D[curFocusWin]) ? 1 : 0;
		rx = (float)cameraInitLocs[sceneIDX].x;
		ry = (float)cameraInitLocs[sceneIDX].y;
		dz = (float)cameraInitLocs[sceneIDX].z;
		setFocus();
	}
	
	//////////////////////////////////////////////////////
	/// user interaction
	//////////////////////////////////////////////////////	
	
	public void keyPressed(){
		switch (key){
			case '1' : {break;}
			case '2' : {break;}
			case '3' : {break;}
			case '4' : {break;}
			case '5' : {break;}							
			case '6' : {break;}
			case '7' : {break;}
			case '8' : {break;}
			case '9' : {break;}
			case '0' : {setFlags(showUIMenu,true); break;}							//to force show UI menu
			case ' ' : {setFlags(runSim,!flags[runSim]); break;}							//run sim
			case 'f' : {setCamView();break;}//reset camera
			case 'a' :
			case 'A' : {setFlags(saveAnim,!flags[saveAnim]);break;}						//start/stop saving every frame for making into animation
			case 's' :
			case 'S' : {save(sketchPath() + "\\"+prjNmLong+dateStr+"\\"+prjNmShrt+"_img"+timeStr + ".jpg");break;}//save picture of current image			
//				case ';' :
//				case ':' : {((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiObjs[((mySideBarMenu)dispWinFrames[dispMenuIDX]).gIDX_cycModDraw].modVal(-1); break;}//decrease the number of cycles between each draw, to some lower bound
//				case '\'' :
//				case '"' : {((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiObjs[((mySideBarMenu)dispWinFrames[dispMenuIDX]).gIDX_cycModDraw].modVal(1); break;}//increase the number of cycles between each draw to some upper bound		
			default : {	}
		}//switch	
		
		if((!flags[shiftKeyPressed])&&(key==CODED)){setFlags(shiftKeyPressed,(keyCode  == KeyEvent.VK_SHIFT));}
		if((!flags[altKeyPressed])&&(key==CODED)){setFlags(altKeyPressed,(keyCode  == KeyEvent.VK_ALT));}
		if((!flags[cntlKeyPressed])&&(key==CODED)){setFlags(cntlKeyPressed,(keyCode  == KeyEvent.VK_CONTROL));}
	}
	public void keyReleased(){
		if((flags[shiftKeyPressed])&&(key==CODED)){ if(keyCode == KeyEvent.VK_SHIFT){endShiftKey();}}
		if((flags[altKeyPressed])&&(key==CODED)){ if(keyCode == KeyEvent.VK_ALT){endAltKey();}}
		if((flags[cntlKeyPressed])&&(key==CODED)){ if(keyCode == KeyEvent.VK_CONTROL){endCntlKey();}}
	}		
	public void endShiftKey(){
		clearFlags(new int []{shiftKeyPressed, modView});
		for(int i =0; i<SnowGlobeWin.numDispWins; ++i){dispWinFrames[i].endShiftKey();}
	}
	public void endAltKey(){
		clearFlags(new int []{altKeyPressed});
		for(int i =0; i<SnowGlobeWin.numDispWins; ++i){dispWinFrames[i].endAltKey();}			
	}
	public void endCntlKey(){
		clearFlags(new int []{cntlKeyPressed});
		for(int i =0; i<SnowGlobeWin.numDispWins; ++i){dispWinFrames[i].endCntlKey();}			
	}

	//2d range checking of point
	public boolean ptInRange(double x, double y, double minX, double minY, double maxX, double maxY){return ((x > minX)&&(x <= maxX)&&(y > minY)&&(y <= maxY));}	
	//gives multiplier based on whether shift, alt or cntl (or any combo) is pressed
	public double clickValModMult(){
		return ((flags[altKeyPressed] ? .1 : 1.0) * (flags[cntlKeyPressed] ? 10.0 : 1.0));			
	}
	
	public void mouseMoved(){for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseMove(mouseX, mouseY,c.getMseLoc(sceneCtrVals[sceneIDX]))){return;}}}
	public void mousePressed() {
		//verify left button if(mouseButton == LEFT)
		setFlags(mouseClicked, true);
		if(mouseButton == LEFT){			mouseClicked(0);} 
		else if (mouseButton == RIGHT) {	mouseClicked(1);}
	}// mousepressed	
	
	private void mouseClicked(int mseBtn){ for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseClick(mouseX, mouseY,c.getMseLoc(sceneCtrVals[sceneIDX]),mseBtn)){return;}}}	
	public void mouseDragged(){//pmouseX is previous mouse x
		if((flags[shiftKeyPressed]) && (canMoveView[curFocusWin])){		//modifying view - always bypass HUD windows if doing this
			flags[modView]=true;
			if(mouseButton == LEFT){			rx-=PI*(mouseY-pmouseY)/height; ry+=PI*(mouseX-pmouseX)/width;} 
			else if (mouseButton == RIGHT) {	dz-=camZScaleMult * (float)(mouseY-pmouseY);}
		} else {
			if(mouseButton == LEFT){			mouseDragged(0);}//mouseLeftDragged();} 
			else if (mouseButton == RIGHT) {	mouseDragged(1);}
			//for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseDrag(mouseX, mouseY, pmouseX, pmouseY,new myVector(c.getOldMseLoc(),c.getMseLoc()))) {return;}}
		}
	}//mouseDragged()
	private void mouseDragged(int mseBtn){
		for(int i =0; i<numDispWins; ++i){if (dispWinFrames[i].handleMouseDrag(mouseX, mouseY, pmouseX, pmouseY,c.getMseLoc(sceneCtrVals[sceneIDX]),new myVector(c.getOldMseLoc(),c.getMseLoc()),mseBtn)) {return;}}		
	}

	public void mouseReleased(){
		clearFlags(new int[]{mouseClicked, modView});
		msClkStr = "";
		for(int i =0; i<numDispWins; ++i){dispWinFrames[i].handleMouseRelease();}
		flags[drawing] = false;
		//c.clearMsDepth();
	}//mouseReleased

	//these tie using the UI buttons to modify the window in with using the boolean tags - PITA but currently necessary
	public void handleShowWin(int btn, int val){handleShowWin(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleShowWin(int btn, int val, boolean callFlags){//{"Score","Curve","InstEdit"},					//display specific windows - multi-select/ always on if sel
		if(!callFlags){//called from setflags - only sets button state in UI
			((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[mySideBarMenu.btnShowWinIdx][btn] = val;
		} else {//called from clicking on buttons in UI
			boolean bVal = (val == 1?  false : true);
			switch(btn){
				case 0 : {setFlags(showSnowGlbWin, bVal);break;}
			}
		}
	}//handleShowWin
	
//		//process request to add a  new component
//		public void handleAddNewCmp(int btn, int val){handleAddNewCmp(btn, val, true);}					//display specific windows - multi-select/ always on if sel
//		public void handleAddNewCmp(int btn, int val, boolean callFlags){//{"Score","Staff","Measure","Note"},			//add - momentary
//			if(!callFlags){
//				((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[mySideBarMenu.btnAddNewCmpIdx][btn] = val;
//			} else {
//				switch(btn){
//					case 0 : {break;}//new song
//					case 1 : {break;}//
//					case 2 : {break;}
//					case 3 : {break;}
//				}
//			}
//		}//handleAddNewCmp
	
	//process to delete an existing component
	public void handleDBGSelCmp(int btn, int val){handleDBGSelCmp(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleDBGSelCmp(int btn, int val, boolean callFlags){//{"Score","Staff","Measure","Note"},			//del - momentary
		if(!callFlags){
			((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[mySideBarMenu.btnDBGSelCmpIdx][btn] = val;
		} else {
			switch(btn){
				case 0 : {dispWinFrames[curFocusWin].clickDebug(btn) ;break;}
				case 1 : {dispWinFrames[curFocusWin].clickDebug(btn) ;break;}
				case 2 : {dispWinFrames[curFocusWin].clickDebug(btn) ;break;}
				case 3 : {dispWinFrames[curFocusWin].clickDebug(btn) ;break;}
			}
		}
	}//handleAddDelSelCmp	
	
//		//process to edit an instrument : TODO
//		public void handleInstEdit(int btn, int val){handleInstEdit(btn, val, true);}					//display specific windows - multi-select/ always on if sel
//		public void handleInstEdit(int btn, int val, boolean callFlags){//{"New","Edit","Delete"},					//Instrument edit/modify - momentary
//			if(!callFlags){
//				((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[mySideBarMenu.btnInstEditIdx][btn] = val;
//			} else {
//				switch(btn){
//					case 0 : {break;}
//					case 1 : {break;}
//					case 2 : {break;}
//				}		
//			}
//		}//handleInstEdit		
//		
	//process to handle file io		
	public void handleFileCmd(int btn, int val){handleFileCmd(btn, val, true);}					//display specific windows - multi-select/ always on if sel
	public void handleFileCmd(int btn, int val, boolean callFlags){//{"Load","Save"},							//load an existing score, save an existing score - momentary	
		if(!callFlags){
			((mySideBarMenu)dispWinFrames[dispMenuIDX]).guiBtnSt[mySideBarMenu.btnFileCmdIdx][btn] = val;
		} else {
			switch(btn){
				case 0 : {break;}
				case 1 : {break;}
			}		
		}
	}//handleFileCmd
	public void initDispWins(){
		//float InstEditWinHeight = InstEditWinYMult * height;		//how high is the InstEdit window when shown
		//instanced window dimensions when open and closed - only showing 1 open at a time
		winRectDimOpen[dispSnowGlobeIDX] =  new float[]{menuWidth+hideWinWidth, 0,width-menuWidth-hideWinWidth,height-hidWinHeight};			
		//hidden
		winRectDimClose[dispSnowGlobeIDX] =  new float[]{menuWidth, 0, hideWinWidth, height};				
		
		winTrajFillClrs = new int []{gui_Black,gui_LightGray,gui_LightGray,gui_LightGray,gui_LightGreen};		//set to color constants for each window
		winTrajStrkClrs = new int []{gui_Black,gui_DarkGray,gui_DarkGray,gui_DarkGray,gui_White};		//set to color constants for each window			
		
		String[] winTitles = new String[]{"","Sphere UI","JCuda Sim Window", "Pop-up Edit"},
				winDescr = new String[] {"","Control quantities using the spheres","Simulation Using JCuda", "Trajectory-based Edit Window"};
//			//display window initialization	
		int wIdx = dispSnowGlobeIDX , fIdx = showSnowGlbWin;
		dispWinFrames[wIdx] = new mySnowGlobeWin(this, winTitles[wIdx], fIdx,winFillClrs[wIdx], winStrkClrs[wIdx], winRectDimOpen[wIdx], winRectDimClose[wIdx], winDescr[wIdx],canDrawInWin[wIdx]);
		
		for(int i =0; i < numDispWins; ++i){
			dispWinFrames[i].initDrwnTrajs();
			dispWinFrames[i].dispFlags[myDispWindow.is3DWin] = dispWinIs3D[i];
			dispWinFrames[i].setTrajColors(winTrajFillClrs[i], winTrajStrkClrs[i]);
		}	
		
	}//initDispWins

	
	//get the ui rect values of the "master" ui region (another window) -> this is so ui objects of one window can be made, clicked, and shown displaced from those of the parent windwo
	public double[] getUIRectVals(int idx){
		switch(idx){
		case dispMenuIDX 		: {return new double[0];}			//idx 0 is parent menu sidebar
		case dispSnowGlobeIDX 	: {return dispWinFrames[dispMenuIDX].uiClkCoords;}
		default :  return dispWinFrames[dispMenuIDX].uiClkCoords;
		}
	}//getUIRectVals	
	
	
//////////////////////////////////////////
/// graphics and base functionality utilities and variables
//////////////////////////////////////////
	//constant path strings for different file types
	public static final String fileDelim = "\\";	
	//display-related size variables
	public final int grid2D_X=800, grid2D_Y=800;	
	public final int gridDimX = 800, gridDimY = 800, gridDimZ = 800;				//dimensions of 3d region
	
	public int scrWidth, scrHeight;			//set to be applet.width and applet.height unless otherwise specified below
	public final int scrWidthMod = 200, 
			scrHeightMod = 0;
	public final float frate = 120;			//frame rate - # of playback updates per second
	
	public int sceneIDX;			//which idx in the 2d arrays of focus vals and glbl center vals to use, based on state
	public myVector[] sceneFcsVals = new myVector[]{						//set these values to be different targets of focus
			new myVector(-grid2D_X/2,-grid2D_Y/1.75f,0),
			new myVector(0,0,0)
	};
	
	public myPoint[] sceneCtrVals = new myPoint[]{						//set these values to be different display center translations -
			new myPoint(0,0,0),										// to be used to calculate mouse offset in world for pick
			new myPoint(0,0,0),										// to be used to calculate mouse offset in world for pick
			new myPoint(-gridDimX/2.0,-gridDimY/2.0,-gridDimZ/2.0)
	};
	
	public float[] camVals, perspVals;		

	public final float 
				camZScaleMult = 2.0f,
				camInitialDist = -.5f*groundRadius*camZScaleMult,		//initial distance camera is from scene - needs to be negative
				camInitRy = 0,
				camInitRx = -PI/2.0f;
	private float dz=0, rx=-0.06f*TWO_PI, ry=-0.04f*TWO_PI;		// distance to camera. Manipulated with wheel or when,view angles manipulated when space pressed but not mouse	
	
	public myVector[] cameraInitLocs = new myVector[]{						//set these values to be different initial camera locations based on 2d or 3d
			new myVector(camInitRx,camInitRy,camInitialDist),
			new myVector(camInitRx*.25,camInitRy,camInitialDist),//<--globe uses this one 
			new myVector(-0.47f,-0.61f,-gridDimZ*2.5f)			
		};
	
	//static variables - put obj constructor counters here
	public static int GUIObjID = 0;										//counter variable for gui objs
	
	//visualization variables
	// boolean flags used to control various elements of the program 
	public boolean[] flags;
	//dev/debug flags
	public final int debugMode 			= 0;			//whether we are in debug mode or not	
	public final int saveAnim 			= 1;			//whether we are saving or not
	//interface flags	
	public final int shiftKeyPressed 	= 2;			//shift pressed
	public final int altKeyPressed  	= 3;			//alt pressed
	public final int cntlKeyPressed  	= 4;			//cntrl pressed
	public final int mouseClicked 		= 5;			//mouse left button is held down	
	public final int drawing			= 6; 			//currently drawing
	public final int modView	 		= 7;			//shift+mouse click+mouse move being used to modify the view
	
	public final int runSim				= 8;			//run simulation (if off localization progresses on single pose		
	public final int showUIMenu 		= 9;			//whether or not to show sidebar menu

	public final int showSnowGlbWin		= 10;
	//public final int showCudaWin		= 11;
	//public final int showPopEdWin		= 12;			
	
	public final int flipDrawnTraj  	= 11;			//whether or not to flip the direction of the drawn melody trajectory
	public final int clearTrajWithNew	= 12;			//clear out existing trajectories and their data when new ones are drawn
	
	public final int numFlags = 13;

	//flags to actually display in menu as clickable text labels - order does matter
	public List<Integer> flagsToShow = Arrays.asList( 
			debugMode, 			
			saveAnim,	
			runSim,
			clearTrajWithNew
			);
	
	public final int numFlagsToShow = flagsToShow.size();
	
	public List<Integer> stateFlagsToShow = Arrays.asList( 
			 shiftKeyPressed,			//shift pressed
			 altKeyPressed  ,			//alt pressed
			 cntlKeyPressed ,			//cntrl pressed
			 mouseClicked 	,			//mouse left button is held down	
			 drawing		, 			//currently drawing
			 modView	 				//shift+mouse click+mouse move being used to modify the view					
			);
	public final int numStFlagsToShow = stateFlagsToShow.size();	

	
	//individual display/HUD windows for gui/user interaction
	public myDispWindow[] dispWinFrames;
	//idx's in dispWinFrames for each window
	public static final int dispMenuIDX = 0,
							dispSnowGlobeIDX = 1;
	
	public static final int numDispWins = 2;	
			
	public int curFocusWin;				//which myDispWindow currently has focus 
	
	//whether or not the display windows will accept a drawn trajectory -- eventually set InstEdit window to be drawable
	public boolean[] canDrawInWin = new boolean[]{	false,true};		
	public boolean[] canShow3DBox = new boolean[]{	false,false};		
	public boolean[] canMoveView = new boolean[]{	false,true};		
	public static final boolean[] dispWinIs3D = new boolean[]{false,true};
	
	public static final int[][] winFillClrs = new int[][]{          
		new int[]{255,255,255,255},                                 	// dispMenuIDX = 0,
		new int[]{0,0,0,255}                                        	// dispSnoGlobe = 2;
	};
	public static final int[][] winStrkClrs = new int[][]{
		new int[]{0,0,0,255},                                    		//dispMenuIDX = 0,
		new int[]{255,255,255,255}                               		//dispSnoGlobe = 2;
	};
	public static int[] winTrajFillClrs = new int []{0,0};		//set to color constants for each window
	public static int[] winTrajStrkClrs = new int []{0,0};		//set to color constants for each window

	
	//unblocked window dimensions - location and dim of window if window is one\
	public float[][] winRectDimOpen;// = new float[][]{new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0}};
	//window dimensions if closed -location and dim of all windows if this window is closed
	public float[][] winRectDimClose;// = new float[][]{new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0},new float[]{0,0,0,0}};
	
	public boolean showInfo;										//whether or not to show start up instructions for code		
	public myVector focusTar;										//target of focus - used in translate to set where the camera is looking - 
																	//set array of vector values (sceneFcsVals) based on application
	//private boolean cyclModCmp;										//comparison every draw of cycleModDraw			
	public final int[] bground = new int[]{244,244,244,255};		//bground color
			
	public myPoint mseCurLoc2D;
	//how many frames to wait to actually refresh/draw
	//public int cycleModDraw = 1;
	public final int maxCycModDraw = 20;	//max val for cyc mod draw		

	// path and filename to save pictures for animation
	public String animPath, animFileName;
	public int animCounter;	
	public final int scrMsgTime = 50;									//5 seconds to delay a message 60 fps (used against draw count)
	public ArrayDeque<String> consoleStrings;							//data being printed to console - show on screen
	
	public int drawCount,simCycles;												// counter for draw cycles		
	public float menuWidth,menuWidthMult = .15f, hideWinWidth, hideWinWidthMult = .03f, hidWinHeight, hideWinHeightMult = .05f;			//side menu is 15% of screen grid2D_X, 

	public ArrayList<String> DebugInfoAra;										//enable drawing dbug info onto screen
	public String debugInfoString;
	
	//animation control variables	
	public float animCntr = 0, animModMult = 1.0f;
	public final float maxAnimCntr = PI*1000.0f, baseAnimSpd = 1.0f;

	my3DCanvas c;												//3d interaction stuff and mouse tracking
	public String dateStr, timeStr;								//used to build directory and file names for screencaps
	
	public PGraphicsOpenGL pg; 
	public PGL pgl;
	//public GL2 gl;
	
	public final double epsVal = .000000001, msClkEps = 40;				//calc epsilon, distance within which to check if clicked from a point
	public float feps = .000001f;
	public float SQRT2 = sqrt(2.0f);

	public int[] rgbClrs = new int[]{gui_Red,gui_Green,gui_Blue};
	//3dbox stuff
	public myVector[] boxNorms = new myVector[] {new myVector(1,0,0),new myVector(-1,0,0),new myVector(0,1,0),new myVector(0,-1,0),new myVector(0,0,1),new myVector(0,0,-1)};//normals to 3 d bounding boxes
	private final float hGDimX = gridDimX/2.0f, hGDimY = gridDimY/2.0f, hGDimZ = gridDimZ/2.0f;
	private final float tGDimX = gridDimX*10, tGDimY = gridDimY*10, tGDimZ = gridDimZ*20;
	public myPoint[][] boxWallPts = new myPoint[][] {//pts to check if intersection with 3D bounding box happens
			new myPoint[] {new myPoint(hGDimX,tGDimY,tGDimZ), new myPoint(hGDimX,-tGDimY,tGDimZ), new myPoint(hGDimX,tGDimY,-tGDimZ)  },
			new myPoint[] {new myPoint(-hGDimX,tGDimY,tGDimZ), new myPoint(-hGDimX,-tGDimY,tGDimZ), new myPoint(-hGDimX,tGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,hGDimY,tGDimZ), new myPoint(-tGDimX,hGDimY,tGDimZ), new myPoint(tGDimX,hGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,-hGDimY,tGDimZ),new myPoint(-tGDimX,-hGDimY,tGDimZ),new myPoint(tGDimX,-hGDimY,-tGDimZ) },
			new myPoint[] {new myPoint(tGDimX,tGDimY,hGDimZ), new myPoint(-tGDimX,tGDimY,hGDimZ), new myPoint(tGDimX,-tGDimY,hGDimZ)  },
			new myPoint[] {new myPoint(tGDimX,tGDimY,-hGDimZ),new myPoint(-tGDimX,tGDimY,-hGDimZ),new myPoint(tGDimX,-tGDimY,-hGDimZ)  }
	};
	//for multithreading 
	public ExecutorService th_exec;
	public int numThreadsAvail;	

	
	//precalculated values
	public static float piO2 = (float) (Math.PI/2.0), piO7 = (float) (Math.PI/7.0), sPiO7 = (float) (Math.sin(piO7)),cPiO7 = (float) (Math.cos(piO7)), 
						piO6 = (float) (Math.PI/6.0), sPiO6 = (float) (Math.sin(piO6)), cPiO6 = (float) (Math.cos(piO6)), scPiO6 = (float) (sPiO6 * cPiO6),
						piO12 = (float) (Math.PI/12.0), sPiO12 = (float) (Math.sin(piO12)), cPiO12 = (float) (Math.cos(piO12));
	//precalced sin and cos
	//precalc cos and sin
	public float[] cosAra,sinAra;
	public final int numThVals = 43008;//even multiples of potential partitions of sin or cos ara
	
	
///////////////////////////////////
/// generic graphics functions and classes
///////////////////////////////////
		//1 time initialization of things that won't change
	public void initVisOnce(){	
		initProcVals();
		dateStr = "_"+day() + "-"+ month()+ "-"+year();
		timeStr = "_"+hour()+"-"+minute()+"-"+second();
		
		scrWidth = width + scrWidthMod;
		scrHeight = height + scrHeightMod;		//set to be applet.width and applet.height unless otherwise specified below
		
		consoleStrings = new ArrayDeque<String>();				//data being printed to console		
		menuWidth = width * menuWidthMult;						//grid2D_X of menu region	
		hideWinWidth = width * hideWinWidthMult;				//dims for hidden windows
		hidWinHeight = height * hideWinHeightMult;
		c = new my3DCanvas(this);			
		winRectDimOpen = new float[numDispWins][];
		winRectDimClose = new float[numDispWins][];
		winRectDimOpen[0] =  new float[]{0,0, menuWidth, height};
		winRectDimClose[0] =  new float[]{0,0, hideWinWidth, height};
		
		//display window initialization
		dispWinFrames = new myDispWindow[numDispWins];		
		//menu bar init
		dispWinFrames[dispMenuIDX] = new mySideBarMenu(this, "UI Window", showUIMenu,  winFillClrs[dispMenuIDX], winStrkClrs[dispMenuIDX], winRectDimOpen[dispMenuIDX],winRectDimClose[dispMenuIDX], "User Controls",canDrawInWin[dispMenuIDX]);			
		
		colorMode(RGB, 255, 255, 255, 255);
		mseCurLoc2D = new myPoint(0,0,0);	
		initBoolFlags();
		float camZ = ((height/2.0f) / tan(PI/6.0f))*this.camZScaleMult;
		camVals = new float[]{width/2.0f, height/2.0f, camZ, width/2.0f, height/2.0f, 0, 0, 1, 0};
		perspVals = new float[]{PI/3.0f, width/height, camZ*.1f, camZ*10.0f};
		showInfo = true;
		outStr2Scr("Current sketchPath " + sketchPath());
		
		initCamView();
		simCycles = 0;
		animPath = sketchPath() + "\\"+prjNmLong+"_" + (int) random(1000);
		animFileName = "\\" + prjNmLong;
		
		//precalc cos and sin
		cosAra = new float[numThVals];
		sinAra = new float[numThVals];
		for (int i = 0; i < numThVals; ++i) {
			float theta = (TWO_PI * i )/numThVals;
			cosAra[i] = cos(theta);
			sinAra[i] = sin(theta);
		}
	}				
	
	//initialize processing values
	private void initProcVals(){
		frameRate(frate);		
		strokeCap(SQUARE);//makes the ends of stroke lines squared off		
		sphereDetail(4);	
		textureMode(NORMAL);			
		rectMode(CORNER);			
	}
	
		//init boolean state machine flags for program
	public void initBoolFlags(){
		flags = new boolean[numFlags];
		for (int i = 0; i < numFlags; ++i) { flags[i] = false;}	
		((mySideBarMenu)dispWinFrames[dispMenuIDX]).initPFlagColors();			//init sidebar window flags
	}		
	
	//address all flag-setting here, so that if any special cases need to be addressed they can be
	public void setFlags(int idx, boolean val ){
		flags[idx] = val;
		switch (idx){
			case debugMode 			: { break;}//anything special for debugMode 			
			case saveAnim 			: { break;}//anything special for saveAnim 			
			case altKeyPressed 		: { break;}//anything special for altKeyPressed 	
			case shiftKeyPressed 	: { break;}//anything special for shiftKeyPressed 	
			case mouseClicked 		: { break;}//anything special for mouseClicked 		
			case modView	 		: { break;}//anything special for modView	 	
			case drawing			: { break;}
			case runSim			: {break;}// handleTrnsprt((val ? 2 : 1) ,(val ? 1 : 0),false); break;}		//anything special for runSim	
			//case flipDrawnTraj		: { dispWinFrames[dispPianoRollIDX].rebuildDrawnTraj();break;}						//whether or not to flip the drawn melody trajectory, width-wise
			case flipDrawnTraj		: { for(int i =1; i<dispWinFrames.length;++i){dispWinFrames[i].rebuildAllDrawnTrajs();}break;}						//whether or not to flip the drawn melody trajectory, width-wise
			case showUIMenu 	    : { dispWinFrames[dispMenuIDX].setShow(val);    break;}											//whether or not to show the main ui window (sidebar)
			
			case showSnowGlbWin		: {setWinFlagsXOR(dispSnowGlobeIDX, val); break;}

			//case useDrawnVels 		: {for(int i =1; i<dispWinFrames.length;++i){dispWinFrames[i].rebuildAllDrawnTrajs();}break;}
			default : {break;}
		}
	}//setFlags  
	
	//set the height of each window that is above the InstEdit window, to move up or down when it changes size
	//specify mutually exclusive flags here
	public int[] winFlagsXOR = new int[]{showSnowGlbWin};//showSequence,showSphereUI};
	//specify windows that cannot be shown simultaneously here
	public int[] winDispIdxXOR = new int[]{dispSnowGlobeIDX};//dispPianoRollIDX,dispSphereUIIDX};
	public void setWinFlagsXOR(int idx, boolean val){
		//outStr2Scr("SetWinFlagsXOR : idx " + idx + " val : " + val);
		if(val){//turning one on
			//turn off not shown, turn on shown				
			for(int i =0;i<winDispIdxXOR.length;++i){//skip first window - ui menu - and last window - InstEdit window
				if(winDispIdxXOR[i]!= idx){dispWinFrames[winDispIdxXOR[i]].setShow(false);handleShowWin(i ,0,false); flags[winFlagsXOR[i]] = false;}
				else {
					dispWinFrames[idx].setShow(true);
					handleShowWin(i ,1,false); 
					flags[winFlagsXOR[i]] = true;
					curFocusWin = winDispIdxXOR[i];
					setCamView();
				}
			}
		} else {				//if turning off a window - need a default uncloseable window - for now just turn on next window : idx-1 is idx of allowable winwdows (idx 0 is sidebar menu)
			setWinFlagsXOR((((idx-1) + 1) % winFlagsXOR.length)+1, true);
		}			
	}//setWinFlagsXOR
	
	//set flags appropriately when only 1 can be true 
	public void setFlagsXOR(int tIdx, int[] fIdx){for(int i =0;i<fIdx.length;++i){if(tIdx != fIdx[i]){flags[fIdx[i]] =false;}}}				
	public void clearFlags(int[] idxs){		for(int idx : idxs){flags[idx]=false;}	}			
		//called every time re-initialized
	public void initVisProg(){	drawCount = 0;		debugInfoString = "";		reInitInfoStr();}	
	public void initCamView(){	dz=camInitialDist;	ry=camInitRy;	rx=camInitRx - ry;	}
	public void reInitInfoStr(){		DebugInfoAra = new ArrayList<String>();		DebugInfoAra.add("");	}	
	public int addInfoStr(String str){return addInfoStr(DebugInfoAra.size(), str);}
	public int addInfoStr(int idx, String str){	
		int lstIdx = DebugInfoAra.size();
		if(idx >= lstIdx){		for(int i = lstIdx; i <= idx; ++i){	DebugInfoAra.add(i,"");	}}
		setInfoStr(idx,str);	return idx;
	}
	public void setInfoStr(int idx, String str){DebugInfoAra.set(idx,str);	}
	public void drawInfoStr(float sc){//draw text on main part of screen
		pushMatrix();		pushStyle();
		fill(0,0,0,100);
		translate((menuWidth),0);
		scale(sc,sc);
		for(int i = 0; i < DebugInfoAra.size(); ++i){		text((flags[debugMode]?(i<10?"0":"")+i+":     " : "") +"     "+DebugInfoAra.get(i)+"\n\n",0,(10+(12*i)));	}
		popStyle();	popMatrix();
	}		
	//vector and point functions to be compatible with earlier code from jarek's class or previous projects	
	//draw bounding box for 3d
	public void drawBoxBnds(){
		pushMatrix();
		pushStyle();
		strokeWeight(3f);
		noFill();
		setColorValStroke(gui_TransGray);
		
		box(gridDimX,gridDimY,gridDimZ);
		popStyle();		
		popMatrix();
	}		
	//drawsInitial setup for each draw
	public void drawSetup(){			
		//perspective(perspVals[0],perspVals[1],perspVals[2],perspVals[3]);
		scale(1.0f/camZScaleMult);
		camera(camVals[0],camVals[1],camVals[2],camVals[3],camVals[4],camVals[5],camVals[6],camVals[7],camVals[8]);       
		//outStr2Scr("rx :  " + rx + " ry : " + ry + " dz : " + dz);
		translate((float)width/2.0f,(float)height/2.0f,(float)dz); // puts origin of model at screen center and moves forward/away by dz
		scale(camZScaleMult);
	    setCamOrient();
	    turnOnLights();
	}//drawSetup	
	//turn on lights for this sketch
	public void turnOnLights(){
	    lights(); 
	}
	public void setCamOrient(){rotateX(rx);rotateY(ry); rotateX(PI/(2.0f));		}//sets the rx, ry, pi/2 orientation of the camera eye	
	public void unSetCamOrient(){rotateX(-PI/(2.0f)); rotateY(-ry);   rotateX(-rx); }//reverses the rx,ry,pi/2 orientation of the camera eye - paints on screen and is unaffected by camera movement
	public void drawAxes(double len, float stW, myPoint ctr, int alpha, boolean centered){//axes using current global orientation
		pushMatrix();pushStyle();
			strokeWeight(stW/camZScaleMult);
			stroke(255,0,0,alpha);
			if(centered){
				double off = len*.5f;
				line(ctr.x-off,ctr.y,ctr.z,ctr.x+off,ctr.y,ctr.z);stroke(0,255,0,alpha);line(ctr.x,ctr.y-off,ctr.z,ctr.x,ctr.y+off,ctr.z);stroke(0,0,255,alpha);line(ctr.x,ctr.y,ctr.z-off,ctr.x,ctr.y,ctr.z+off);} 
			else {		line(ctr.x,ctr.y,ctr.z,ctr.x+len,ctr.y,ctr.z);stroke(0,255,0,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x,ctr.y+len,ctr.z);stroke(0,0,255,alpha);line(ctr.x,ctr.y,ctr.z,ctr.x,ctr.y,ctr.z+len);}
		popStyle();	popMatrix();	
	}//	drawAxes
	public void drawAxes(double len, float stW, myPoint ctr, myVector[] _axis, int alpha, boolean drawVerts){//RGB -> XYZ axes
		pushMatrix();pushStyle();
		if(drawVerts){
			show(ctr,3,gui_Black,gui_Black, false);
			for(int i=0;i<_axis.length;++i){show(myPoint._add(ctr, myVector._mult(_axis[i],len)),3,rgbClrs[i],rgbClrs[i], false);}
		}
		strokeWeight(stW/camZScaleMult);
		for(int i =0; i<3;++i){	setColorValStroke(rgbClrs[i]);	showVec(ctr,len, _axis[i]);	}
		popStyle();	popMatrix();	
	}//	drawAxes
	public void drawAxes(double len, float stW, myPoint ctr, myVector[] _axis, int[] clr, boolean drawVerts){//all axes same color
		pushMatrix();pushStyle();
			if(drawVerts){
				show(ctr,2,gui_Black,gui_Black, false);
				for(int i=0;i<_axis.length;++i){show(myPoint._add(ctr, myVector._mult(_axis[i],len)),2,rgbClrs[i],rgbClrs[i], false);}
			}
			strokeWeight(stW/camZScaleMult);stroke(clr[0],clr[1],clr[2],clr[3]);
			for(int i =0; i<3;++i){	showVec(ctr,len, _axis[i]);	}
		popStyle();	popMatrix();	
	}//	drawAxes

	public void drawText(String str, double x, double y, double z, int clr){
		int[] c = getClr(clr);
		pushMatrix();	pushStyle();
			fill(c[0],c[1],c[2],c[3]);
			unSetCamOrient();
			translate((float)x,(float)y,(float)z);
			text(str,0,0,0);		
		popStyle();	popMatrix();	
	}//drawText	
	public void savePic(){		save(animPath + animFileName + ((animCounter < 10) ? "000" : ((animCounter < 100) ? "00" : ((animCounter < 1000) ? "0" : ""))) + animCounter + ".jpg");		animCounter++;		}
	public void line(double x1, double y1, double z1, double x2, double y2, double z2){line((float)x1,(float)y1,(float)z1,(float)x2,(float)y2,(float)z2 );}
	public void line(myPoint p1, myPoint p2){line((float)p1.x,(float)p1.y,(float)p1.z,(float)p2.x,(float)p2.y,(float)p2.z);}
	public void line(myPointf p1, myPointf p2){line(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);}
	
	//very fast mechanism for setting an array of doubles to passed value - takes advantage of arraycopy copying blocks at a time of a size > 1
	public void dAraFill(double[] ara, double val){
		  int len = ara.length;
		  if (len > 0){ara[0] = val; }
		  for (int i = 1; i < len; i += i){
			  System.arraycopy(ara, 0, ara, i, ((len - i) < i) ? (len - i) : i);
		  }		
	}
	
	public double min3(double x, double y, double z) { return (x < y ? (x < z ? x : z) : (y < z ? y : z)); }
	public double max3(double x, double y, double z) { return (x > y ? (x > z ? x : z) : (y > z ? y : z)); }

	//////////////////
	///pshape cylinder stuff
	//////////////////
	//draws a quad with normals at each pair of verts - used to build smooth cylinder
	public void drawShape(float[] x, float[] y,float[] z){
		beginShape();
		for(int i=0;i<x.length;++i){		
			vertex (x[i],y[i],z[i]);				
		}	
		endShape(CLOSE);	
	}//drawShape
	
	//draws a quad with normals at each pair of verts - used to build smooth cylinder
	public void drawShape(myVectorf[] verts){
		beginShape();
		noStroke();
		//setShClr(sh,setClr);
		for(int i=0;i<verts.length;++i){
			vertex (verts[i].x,verts[i].y,verts[i].z);				
		}
		endShape(CLOSE);
	}//drawShape

	//draws a quad with normals at each pair of verts - used to build smooth cylinder
	public PShape buildVShape(myVectorf[] verts, int setClr){
		PShape sh = createShape();
		sh.beginShape();
		sh.noStroke();
		setShClr(sh,setClr);
		for(int i=0;i<verts.length;++i){
			sh.vertex (verts[i].x,verts[i].y,verts[i].z);				
		}
		sh.endShape(CLOSE);	
		return sh;
	}//drawShape
	
	/**
	*  cylinder with different radii for top and bottom - oriented along y axis
	*/		
	public PShape buildCylinder (float radBot, float radTop, float smHeight, int sides, boolean drawBottom, boolean drawTop, int setClr, myVectorf trans) {
		PShape shRes = createShape(GROUP);
		int ii;
		//how many to increment the cosine and sine arrays
		int incr = numThVals/sides;
		if (drawBottom){shRes.addChild(buildCap(incr, 0, -1.0f,radBot, setClr, trans));}//drawbottom
		if (drawTop){shRes.addChild(buildCap(incr, smHeight, 1.0f,radTop, setClr, trans));}//drawTop	  
		// main body of cylinder
		myVectorf[] verts = new myVectorf[]{new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0)};
		float smHsq = smHeight*smHeight, radDiff = (radTop - radBot),radDiffSq = radDiff*radDiff, hyp = sqrt(smHsq + radDiffSq);
		
		for (int i = 0; i < numThVals; i+=incr) {
			ii = (i+incr) % numThVals;
			verts[0].set((cosAra[i]  * radBot), 0, (sinAra[i] * radBot));
			verts[1].set((cosAra[i]  * radTop), smHeight, (sinAra[i] * radTop));
			verts[2].set((cosAra[ii] * radTop), smHeight, (sinAra[ii] * radTop));
			verts[3].set((cosAra[ii] * radBot), 0, (sinAra[ii] * radBot));
			for(int j =0; j<4;++j){verts[j].set(myVectorf._add(verts[j], trans));}
			shRes.addChild(buildVShape(verts, setClr));
		}
		return shRes;
	}//buildCylinder	
	
	
	private PShape buildCap(int incr, float height, float normDir, float rad, int setClr, myVectorf trans){
		int ii;
		PShape sh = createShape();
		sh.beginShape( TRIANGLES);
		sh.noStroke();
		sh.translate(trans.x, trans.y, trans.z);
		setShClr(sh,setClr);
		if(normDir > 0){
			for (int i = 0; i < numThVals; i+=incr) {
				ii = (i+incr) % numThVals;
				sh.vertex ( (cosAra[ii] * rad), height, (sinAra[ii] * rad));	 sh.vertex ( (cosAra[i] * rad),height, (sinAra[i] * rad));  sh.vertex (0.0f, height, 0.0f);
			}
		}else{
			for (int i = 0; i < numThVals; i+=incr) {
				ii = (i+incr) % numThVals;
				sh.vertex (0.0f, height, 0.0f);sh.vertex ( (cosAra[i] * rad),height, (sinAra[i] * rad));sh.vertex ( (cosAra[ii] * rad), height, (sinAra[ii] * rad));	 
			}
		}
		sh.endShape(PConstants.CLOSE);
		return sh;
	}//drawCap
	
	/**
	*  cylinder with different radii for top and bottom
	*/		
	public void drawCylinder (float radBot, float radTop, float smHeight, int sides, boolean drawBottom, boolean drawTop) {
		pushMatrix();pushStyle();
		int ii;
		//how many to increment the cosine and sine arrays
		int incr = numThVals/sides;
		if (drawBottom){drawCap(incr, 0, -1.0f,radBot);}//drawbottom
		if (drawTop){drawCap(incr, smHeight, 1.0f,radTop);}//drawTop	  
		// main body of cylinder
		myVectorf[] verts = new myVectorf[]{new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0), new myVectorf(0,0,0)};
		for (int i = 0; i < numThVals; i+=incr) {
			ii = (i+incr) % numThVals;
			verts[0].set((cosAra[i]  * radBot), 0, (sinAra[i] * radBot));
			verts[1].set((cosAra[i]  * radTop), smHeight, (sinAra[i] * radTop));
			verts[2].set((cosAra[ii] * radTop), smHeight, (sinAra[ii] * radTop));
			verts[3].set((cosAra[ii] * radBot), 0, (sinAra[ii] * radBot));
		//	for(int j =0; j<4;++j){verts[j].set(myVectorf._add(verts[j], trans));}
			drawShape(verts);
		}

		popStyle();popMatrix();
	}//drawcylinder
	
	private void drawCap(int incr, float height, float normDir, float rad){
		int ii;
		beginShape( PConstants.TRIANGLES);
		if(normDir > 0){
			for (int i = 0; i < numThVals; i+=incr) {
				ii = (i+incr) % numThVals;
				vertex ( (cosAra[ii] * rad), height, (sinAra[ii] * rad));	 vertex ( (cosAra[i] * rad),height, (sinAra[i] * rad));  vertex (0.0f, height, 0.0f);
			}
		}else{
			for (int i = 0; i < numThVals; i+=incr) {
				ii = (i+incr) % numThVals;
				vertex (0.0f, height, 0.0f);vertex ( (cosAra[i] * rad),height, (sinAra[i] * rad));vertex ( (cosAra[ii] * rad), height, (sinAra[ii] * rad));	 
			}
		}
		endShape(PConstants.CLOSE);
	}//drawCap
	
//////////////////
///end cylinder stuff
//////////////////
	
	public void drawOnScreenData(){
		if(flags[debugMode]){
			pushMatrix();pushStyle();			
			reInitInfoStr();
			addInfoStr(0,"mse loc on screen : " + new myPoint(mouseX, mouseY,0) + " mse loc in world :"+c.mseLoc +"  Eye loc in world :"+ c.eyeInWorld); 
			String[] res = ((mySideBarMenu)dispWinFrames[dispMenuIDX]).getDebugData();		//get debug data for each UI object
			//for(int s=0;s<res.length;++s) {	addInfoStr(res[s]);}				//add info to string to be displayed for debug
			int numToPrint = min(res.length,80);
			for(int s=0;s<numToPrint;++s) {	addInfoStr(res[s]);}				//add info to string to be displayed for debug
			drawInfoStr(1.0f); 	
			popStyle();	popMatrix();		
		}
		else if(showInfo){
			pushMatrix();pushStyle();			
			reInitInfoStr();	
			if(showInfo){
//			      addInfoStr(0,"Click the light green box to the left to toggle showing this message.");
//			      addInfoStr(1,"--Shift-Click-Drag to change view.  Shift-RClick-Drag to zoom.");
             // addInfoStr(3,"Values at Mouse Location : "+ values at mouse location);
			}
			String[] res = consoleStrings.toArray(new String[0]);
			int dispNum = min(res.length, 80);
			for(int i=0;i<dispNum;++i){addInfoStr(res[i]);}
		    drawInfoStr(1.1f); 
			popStyle();	popMatrix();	
		}
	}
	//print out multiple-line text to screen
	public void ml_text(String str, float x, float y){
		String[] res = str.split("\\r?\\n");
		float disp = 0;
		for(int i =0; i<res.length; ++i){
			text(res[i],x, y+disp);		//add console string output to screen display- decays over time
			disp += 12;
		}
	}
	
	//print out string in display window
	public void outStr2Scr(String str){outStr2Scr(str,true);}
	//print informational string data to console, and to screen
	public void outStr2Scr(String str, boolean showDraw){
		if(trim(str) != ""){	System.out.println(str);}
		String[] res = str.split("\\r?\\n");
		if(showDraw){
			for(int i =0; i<res.length; ++i){
				consoleStrings.add(res[i]);		//add console string output to screen display- decays over time
			}
		}
	}
	//utilities

	//handle user-driven file load or save - returns a filename + filepath string
	public String FileSelected(File selection){
		if (null==selection){return null;}
		return selection.getAbsolutePath();		
	}//FileSelected

//		//s-cut to print to console
	public void pr(String str){outStr2Scr(str);}
	
	public String getFName(String fNameAndPath){
		String[] strs = fNameAndPath.split("/");
		return strs[strs.length-1];
	}
	
	//load a file as text strings
	public String[] loadFileIntoStringAra(String fileName, String dispYesStr, String dispNoStr){
		String[] strs = null;
		try{
			strs = loadStrings(fileName);
			System.out.println(dispYesStr+"\tLength : " + strs.length);
		} catch (Exception e){System.out.println("!!"+dispNoStr);return null;}
		return strs;		
	}//loadFileIntoStrings

	//public void scribeHeaderRight(String s) {scribeHeaderRight(s, 20);} // writes black on screen top, right-aligned
	//public void scribeHeaderRight(String s, float y) {fill(0); text(s,width-6*s.length(),y); noFill();} // writes black on screen top, right-aligned
	public void displayHeader() { // Displays title and authors face on screen
		hint(PConstants.DISABLE_DEPTH_TEST);
	    float stVal = 17;
	    int idx = 1;	
	    translate(0,10);
	    fill(0); text("Shift-Click-Drag to change view.",width-190, stVal*idx++); noFill(); 
	    fill(0); text("Shift-RClick-Drag to zoom.",width-160, stVal*idx++); noFill();
	    fill(0); text("John Turner",width-75, stVal*idx++); noFill();			
	    //image(jtFace,  width-111,stVal*idx,100,100);
		hint(PConstants.ENABLE_DEPTH_TEST);
	}
	
	//project passed point onto box surface based on location - to help visualize the location in 3d
	public void drawProjOnBox(myPoint p){
		//myPoint[]  projOnPlanes = new myPoint[6];
		myPoint prjOnPlane;
		//public myPoint intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		pushMatrix();
		translate(-p.x,-p.y,-p.z);
		for(int i  = 0; i< 6; ++i){				
			prjOnPlane = bndChkInCntrdBox3D(intersectPl(p, boxNorms[i], boxWallPts[i][0],boxWallPts[i][1],boxWallPts[i][2]));				
			show(prjOnPlane,5,rgbClrs[i/2],rgbClrs[i/2], false);				
		}
		popMatrix();
	}//drawProjOnBox
	
	public myPoint getScrLocOf3dWrldPt(myPoint pt){	return new myPoint(screenX((float)pt.x,(float)pt.y,(float)pt.z),screenY((float)pt.x,(float)pt.y,(float)pt.z),screenZ((float)pt.x,(float)pt.y,(float)pt.z));}
	
	public myPoint bndChkInBox2D(myPoint p){p.set(Math.max(0,Math.min(p.x,grid2D_X)),Math.max(0,Math.min(p.y,grid2D_Y)),0);return p;}
	public myPoint bndChkInBox3D(myPoint p){p.set(Math.max(0,Math.min(p.x,gridDimX)), Math.max(0,Math.min(p.y,gridDimY)),Math.max(0,Math.min(p.z,gridDimZ)));return p;}	
	public myPoint bndChkInCntrdBox3D(myPoint p){
		p.set(Math.max(-hGDimX,Math.min(p.x,hGDimX)), 
				Math.max(-hGDimY,Math.min(p.y,hGDimY)),
				Math.max(-hGDimZ,Math.min(p.z,hGDimZ)));return p;}	
	 
	public void translate(myPoint p){translate((float)p.x,(float)p.y,(float)p.z);}
	public void translate(myVector p){translate((float)p.x,(float)p.y,(float)p.z);}
	public void translate(double x, double y, double z){translate((float)x,(float)y,(float)z);}
	public void translate(double x, double y){translate((float)x,(float)y);}
	public void rotate(float thet, myPoint axis){rotate(thet, (float)axis.x,(float)axis.y,(float)axis.z);}
	public void rotate(float thet, double x, double y, double z){rotate(thet, (float)x,(float)y,(float)z);}
	//************************************************************************
	//**** SPIRAL
	//************************************************************************
	//3d rotation - rotate P by angle a around point G and axis normal to plane IJ
	public myPoint R(myPoint P, double a, myVector I, myVector J, myPoint G) {
		double x= myVector._dot(new myVector(G,P),U(I)), y=myVector._dot(new myVector(G,P),U(J)); 
		double c=Math.cos(a), s=Math.sin(a); 
		double iXVal = x*c-x-y*s, jYVal= x*s+y*c-y;			
		return myPoint._add(P,iXVal,I,jYVal,J); }; 
		
	public cntlPt R(cntlPt P, double a, myVector I, myVector J, myPoint G) {
		double x= myVector._dot(new myVector(G,P),U(I)), y=myVector._dot(new myVector(G,P),U(J)); 
		double c=Math.cos(a), s=Math.sin(a); 
		double iXVal = x*c-x-y*s, jYVal= x*s+y*c-y;		
		return new cntlPt(this, P(P,iXVal,I,jYVal,J), P.r, P.w); }; 
		
	public myPoint PtOnSpiral(myPoint A, myPoint B, myPoint C, double t) {
		//center is coplanar to A and B, and coplanar to B and C, but not necessarily coplanar to A, B and C
		//so center will be coplanar to mp(A,B) and mp(B,C) - use mpCA midpoint to determine plane mpAB-mpBC plane?
		myPoint mAB = new myPoint(A,.5f, B);
		myPoint mBC = new myPoint(B,.5f, C);
		myPoint mCA = new myPoint(C,.5f, A);
		myVector mI = U(mCA,mAB);
		myVector mTmp = myVector._cross(mI,U(mCA,mBC));
		myVector mJ = U(mTmp._cross(mI));	//I and J are orthonormal
		double a =spiralAngle(A,B,B,C); 
		double s =spiralScale(A,B,B,C);
		
		//myPoint G = spiralCenter(a, s, A, B, mI, mJ); 
		myPoint G = spiralCenter(A, mAB, B, mBC); 
		return new myPoint(G, Math.pow(s,t), R(A,t*a,mI,mJ,G));
	  }
	public double spiralAngle(myPoint A, myPoint B, myPoint C, myPoint D) {return myVector._angleBetween(new myVector(A,B),new myVector(C,D));}
	public double spiralScale(myPoint A, myPoint B, myPoint C, myPoint D) {return myPoint._dist(C,D)/ myPoint._dist(A,B);}
	
	public myPoint R(myPoint Q, myPoint C, myPoint P, myPoint R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
		myVector I0=U(C,P), I1=U(C,R), V=new myVector(C,Q); 
		double c=myPoint._dist(I0,I1), s=Math.sqrt(1.-(c*c)); 
		if(Math.abs(s)<0.00001) return Q;
		myVector J0=V(1./s,I1,-c/s,I0);  
		myVector J1=V(-s,I0,c,J0);  
		double x=V._dot(I0), y=V._dot(J0);  
		return P(Q,x,M(I1,I0),y,M(J1,J0)); 
	} 	
	// spiral given 4 points, AB and CD are edges corresponding through rotation
	public myPoint spiralCenter(myPoint A, myPoint B, myPoint C, myPoint D) {         // new spiral center
		myVector AB=V(A,B), CD=V(C,D), AC=V(A,C);
		double m=CD.magn/AB.magn, n=CD.magn*AB.magn;		
		myVector rotAxis = U(AB._cross(CD));		//expect ab and ac to be coplanar - this is the axis to rotate around to find f
		
		myVector rAB = myVector._rotAroundAxis(AB, rotAxis, PConstants.HALF_PI);
		double c=AB._dot(CD)/n, 
				s=rAB._dot(CD)/n;
		double AB2 = AB._dot(AB), a=AB._dot(AC)/AB2, b=rAB._dot(AC)/AB2;
		double x=(a-m*( a*c+b*s)), y=(b-m*(-a*s+b*c));
		double d=1+m*(m-2*c);  if((c!=1)&&(m!=1)) { x/=d; y/=d; };
		return P(P(A,x,AB),y,rAB);
	  }
	
	
	public void cylinder(myPoint A, myPoint B, float r, int c1, int c2) {
		myPoint P = A;
		myVector V = V(A,B);
		myVector I = c.drawSNorm;//U(Normal(V));
		myVector J = U(N(I,V));
		float da = TWO_PI/36;
		beginShape(QUAD_STRIP);
			for(float a=0; a<=TWO_PI+da; a+=da) {fill(c1); gl_vertex(P(P,r*cos(a),I,r*sin(a),J,0,V)); fill(c2); gl_vertex(P(P,r*cos(a),I,r*sin(a),J,1,V));}
		endShape();
	}
	
	//point functions
	public myPoint P() {return new myPoint(); };                                                                          // point (x,y,z)
	public myPoint P(double x, double y, double z) {return new myPoint(x,y,z); };                                            // point (x,y,z)
	public myPoint P(myPoint A) {return new myPoint(A.x,A.y,A.z); };                                                           // copy of point P
	public myPoint P(myPoint A, double s, myPoint B) {return new myPoint(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
	public myPoint L(myPoint A, double s, myPoint B) {return new myPoint(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
	public myPoint P(myPoint A, myPoint B) {return P((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); }                             // (A+B)/2
	public myPoint P(myPoint A, myPoint B, myPoint C) {return new myPoint((A.x+B.x+C.x)/3.0,(A.y+B.y+C.y)/3.0,(A.z+B.z+C.z)/3.0); };     // (A+B+C)/3
	public myPoint P(myPoint A, myPoint B, myPoint C, myPoint D) {return P(P(A,B),P(C,D)); };                                            // (A+B+C+D)/4
	public myPoint P(double s, myPoint A) {return new myPoint(s*A.x,s*A.y,s*A.z); };                                            // sA
	public myPoint A(myPoint A, myPoint B) {return new myPoint(A.x+B.x,A.y+B.y,A.z+B.z); };                                         // A+B
	public myPoint P(double a, myPoint A, double b, myPoint B) {return A(P(a,A),P(b,B));}                                        // aA+bB 
	public myPoint P(double a, myPoint A, double b, myPoint B, double c, myPoint C) {return A(P(a,A),P(b,B,c,C));}                     // aA+bB+cC 
	public myPoint P(double a, myPoint A, double b, myPoint B, double c, myPoint C, double d, myPoint D){return A(P(a,A,b,B),P(c,C,d,D));}   // aA+bB+cC+dD
	public myPoint P(myPoint P, myVector V) {return new myPoint(P.x + V.x, P.y + V.y, P.z + V.z); }                                 // P+V
	public myPoint P(myPoint P, double s, myVector V) {return new myPoint(P.x+s*V.x,P.y+s*V.y,P.z+s*V.z);}                           // P+sV
	public myPoint P(myPoint O, double x, myVector I, double y, myVector J) {return P(O.x+x*I.x+y*J.x,O.y+x*I.y+y*J.y,O.z+x*I.z+y*J.z);}  // O+xI+yJ
	public myPoint P(myPoint O, double x, myVector I, double y, myVector J, double z, myVector K) {return P(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+kZ
	void makePts(myPoint[] C) {for(int i=0; i<C.length; i++) C[i]=P();}

	//draw a circle - JT
	void circle(myPoint P, float r, myVector I, myVector J, int n) {myPoint[] pts = new myPoint[n];pts[0] = P(P,r,U(I));float a = (2*PI)/(1.0f*n);for(int i=1;i<n;++i){pts[i] = R(pts[i-1],a,J,I,P);}pushMatrix(); pushStyle();noFill(); show(pts);popStyle();popMatrix();}; // render sphere of radius r and center P
	
	void circle(myPoint p, float r){ellipse((float)p.x, (float)p.y, r, r);}
	void circle(float x, float y, float r1, float r2){ellipse(x,y, r1, r2);}
		
	void bezier(myPoint A, myPoint B, myPoint C, myPoint D) {bezier((float)A.x,(float)A.y,(float)A.z,(float)B.x,(float)B.y,(float)B.z,(float)C.x,(float)C.y,(float)C.z,(float)D.x,(float)D.y,(float)D.z);} // draws a cubic Bezier curve with control points A, B, C, D
	void bezier(myPoint [] C) {bezier(C[0],C[1],C[2],C[3]);} // draws a cubic Bezier curve with control points A, B, C, D
	myPoint bezierPoint(myPoint[] C, float t) {return P(bezierPoint((float)C[0].x,(float)C[1].x,(float)C[2].x,(float)C[3].x,(float)t),bezierPoint((float)C[0].y,(float)C[1].y,(float)C[2].y,(float)C[3].y,(float)t),bezierPoint((float)C[0].z,(float)C[1].z,(float)C[2].z,(float)C[3].z,(float)t)); }
	myVector bezierTangent(myPoint[] C, float t) {return V(bezierTangent((float)C[0].x,(float)C[1].x,(float)C[2].x,(float)C[3].x,(float)t),bezierTangent((float)C[0].y,(float)C[1].y,(float)C[2].y,(float)C[3].y,(float)t),bezierTangent((float)C[0].z,(float)C[1].z,(float)C[2].z,(float)C[3].z,(float)t)); }

	
	public myPoint Mouse() {return new myPoint(mouseX, mouseY,0);}                                          			// current mouse location
	public myVector MouseDrag() {return new myVector(mouseX-pmouseX,mouseY-pmouseY,0);};                     			// vector representing recent mouse displacement
	
	//public int color(myPoint p){return color((int)p.x,(int)p.z,(int)p.y);}	//needs to be x,z,y for some reason - to match orientation of color frames in z-up 3d geometry
	public int color(myPoint p){return color((int)p.x,(int)p.y,(int)p.z);}	
	
	// =====  vector functions
	public myVector V() {return new myVector(); };                                                                          // make vector (x,y,z)
	public myVector V(double x, double y, double z) {return new myVector(x,y,z); };                                            // make vector (x,y,z)
	public myVector V(myVector V) {return new myVector(V.x,V.y,V.z); };                                                          // make copy of vector V
	public myVector A(myVector A, myVector B) {return new myVector(A.x+B.x,A.y+B.y,A.z+B.z); };                                       // A+B
	public myVector A(myVector U, float s, myVector V) {return V(U.x+s*V.x,U.y+s*V.y,U.z+s*V.z);};                               // U+sV
	public myVector M(myVector U, myVector V) {return V(U.x-V.x,U.y-V.y,U.z-V.z);};                                              // U-V
	public myVector M(myVector V) {return V(-V.x,-V.y,-V.z);};                                                              // -V
	//public myVector V(myVector A, myVector B) {return new myVector((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); }                      // (A+B)/2
	public myVector V(myVector A, float s, myVector B) {return new myVector(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };      // (1-s)A+sB
	public myVector V(myVector A, myVector B, myVector C) {return new myVector((A.x+B.x+C.x)/3.0,(A.y+B.y+C.y)/3.0,(A.z+B.z+C.z)/3.0); };  // (A+B+C)/3
	public myVector V(myVector A, myVector B, myVector C, myVector D) {return V(V(A,B),V(C,D)); };                                         // (A+B+C+D)/4
	public myVector V(double s, myVector A) {return new myVector(s*A.x,s*A.y,s*A.z); };                                           // sA
	public myVector V(double a, myVector A, double b, myVector B) {return A(V(a,A),V(b,B));}                                       // aA+bB 
	public myVector V(double a, myVector A, double b, myVector B, double c, myVector C) {return A(V(a,A,b,B),V(c,C));}                   // aA+bB+cC
	public myVector V(myPoint P, myPoint Q) {return new myVector(P,Q);};                                          // PQ
	public myVector N(myVector U, myVector V) {return V( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                  // UxV cross product (normal to both)
	public myVector N(myPoint A, myPoint B, myPoint C) {return N(V(A,B),V(A,C)); };                                                   // normal to triangle (A,B,C), not normalized (proportional to area)
	public myVector B(myVector U, myVector V) {return U(N(N(U,V),U)); }        

	public myVectorf Vf() {return new myVectorf(); };                                                                          // make vector (x,y,z)
	public myVectorf Vf(float x, double y, double z) {return new myVectorf(x,y,z); };                                            // make vector (x,y,z)
	public myVectorf Vf(myVectorf V) {return new myVectorf(V.x,V.y,V.z); };                                                          // make copy of vector V
	public myVectorf Af(myVectorf A, myVectorf B) {return new myVectorf(A.x+B.x,A.y+B.y,A.z+B.z); };                                       // A+B
	public myVectorf Af(myVectorf U, float s, myVectorf V) {return Vf(U.x+s*V.x,U.y+s*V.y,U.z+s*V.z);};                               // U+sV
	public myVectorf Mf(myVectorf U, myVectorf V) {return Vf(U.x-V.x,U.y-V.y,U.z-V.z);};                                              // U-V
	public myVectorf Mf(myVectorf V) {return Vf(-V.x,-V.y,-V.z);};                                                              // -V
	//public myVectorf Vf(myVectorf A, myVectorf B) {return new myVectorf((A.x+B.x)/2.0,(A.y+B.y)/2.0,(A.z+B.z)/2.0); }                      // (A+B)/2
	public myVectorf Vf(myVectorf A, float s, myVectorf B) {return new myVectorf(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };      // (1-s)A+sB
	public myVectorf Vf(myVectorf A, myVectorf B, myVectorf C) {return new myVectorf((A.x+B.x+C.x)/3.0f,(A.y+B.y+C.y)/3.0f,(A.z+B.z+C.z)/3.0f); };  // (A+B+C)/3
	public myVectorf Vf(myVectorf A, myVectorf B, myVectorf C, myVectorf D) {return Vf(Vf(A,B),Vf(C,D)); };                                         // (A+B+C+D)/4
	public myVectorf Vf(float s, myVectorf A) {return new myVectorf(s*A.x,s*A.y,s*A.z); };                                           // sA
	public myVectorf Vf(float a, myVectorf A, float b, myVectorf B) {return Af(Vf(a,A),Vf(b,B));}                                       // aA+bB 
	public myVectorf Vf(float a, myVectorf A, float b, myVectorf B, float c, myVectorf C) {return Af(Vf(a,A,b,B),Vf(c,C));}                   // aA+bB+cC
	public myVectorf Vf(myPointf P, myPointf Q) {return new myVectorf(P,Q);};                                          // PQ
	public myVectorf Nf(myVectorf U, myVectorf V) {return Vf( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                  // UxV cross product (normal to both)
	public myVectorf Nf(myPointf A, myPointf B, myPointf C) {return Nf(Vf(A,B),Vf(A,C)); };                                                   // normal to triangle (A,B,C), not normalized (proportional to area)
	public myVectorf Bf(myVectorf U, myVectorf V) {return Uf(Nf(Nf(U,V),U)); }        

	//calculate the normal, tangent, binormal components of passed partVec compared to a passed normal (needs to be normalized)
	public myVector[] getVecFrame(myVector partVec, myVector norm) {
		myVector[] result = new myVector[3];//(2, myVector(0, 0, 0));
		result[0] = myVector._mult(norm,(norm._dot(partVec)));//norm dir
		result[1] = myVector._sub(partVec, result[0]);		//tan dir
		result[2] = myVector._cross(result[0], result[1]);
		return result;
	}
	
	public double d(myVector U, myVector V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
	public double dot(myVector U, myVector V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
	public double det2(myVector U, myVector V) {return -U.y*V.x+U.x*V.y; };                                       		// U|V det product
	public double det3(myVector U, myVector V) {double dist = d(U,V); return Math.sqrt(d(U,U)*d(V,V) - (dist*dist)); };                                // U|V det product
	public double m(myVector U, myVector V, myVector W) {return d(U,N(V,W)); };                                                 // (UxV)*W  mixed product, determinant - measures 6x the volume of the parallelapiped formed by myVectortors
	public double m(myPoint E, myPoint A, myPoint B, myPoint C) {return m(V(E,A),V(E,B),V(E,C));}                                    // det (EA EB EC) is >0 when E sees (A,B,C) clockwise
	public double n2(myVector V) {return (V.x*V.x)+(V.y*V.y)+(V.z*V.z);};                                                   // V*V    norm squared
	public double n(myVector V) {return  Math.sqrt(n2(V));};                                                                // ||V||  norm
	public double d(myPoint P, myPoint Q) {return  myPoint._dist(P, Q); };                            // ||AB|| distance
	public double area(myPoint A, myPoint B, myPoint C) {return n(N(A,B,C))/2; };                                               // area of triangle 
	public double volume(myPoint A, myPoint B, myPoint C, myPoint D) {return m(V(A,B),V(A,C),V(A,D))/6; };                           // volume of tet 
	public boolean parallel (myVector U, myVector V) {return n(N(U,V))<n(U)*n(V)*0.00001; }                              // true if U and V are almost parallel
	public double angle(myPoint A, myPoint B, myPoint C){return angle(V(A,B),V(A,C));}												//angle between AB and AC
	public double angle(myPoint A, myPoint B, myPoint C, myPoint D){return angle(U(A,B),U(C,D));}							//angle between AB and CD
	public double angle(myVector U, myVector V){double angle = Math.atan2(n(N(U,V)),d(U,V)),sign = m(U,V,V(0,0,1));if(sign<0){    angle=-angle;}	return angle;}
	public boolean cw(myVector U, myVector V, myVector W) {return m(U,V,W)>0; };                                               // (UxV)*W>0  U,V,W are clockwise
	public boolean cw(myPoint A, myPoint B, myPoint C, myPoint D) {return volume(A,B,C,D)>0; };                                     // tet is oriented so that A sees B, C, D clockwise 
	public boolean projectsBetween(myPoint P, myPoint A, myPoint B) {return dot(V(A,P),V(A,B))>0 && dot(V(B,P),V(B,A))>0 ; };
	public double distToLine(myPoint P, myPoint A, myPoint B) {double res = det3(U(A,B),V(A,P)); return Double.isNaN(res) ? 0 : res; };		//MAY RETURN NAN IF point P is on line
	public myPoint projectionOnLine(myPoint P, myPoint A, myPoint B) {return P(A,dot(V(A,B),V(A,P))/dot(V(A,B),V(A,B)),V(A,B));}
	public boolean isSame(myPoint A, myPoint B) {return (A.x==B.x)&&(A.y==B.y)&&(A.z==B.z) ;}                                         // A==B
	public boolean isSame(myPoint A, myPoint B, double e) {return ((Math.abs(A.x-B.x)<e)&&(Math.abs(A.y-B.y)<e)&&(Math.abs(A.z-B.z)<e));}                   // ||A-B||<e
	
	public myVector W(double s,myVector V) {return V(s*V.x,s*V.y,s*V.z);}                                                      // sV
	public myVectorf Wf(float s,myVectorf V) {return Vf(s*V.x,s*V.y,s*V.z);}                                                      // sV

	public myVector U(myVector v){myVector u = new myVector(v); return u._normalize(); }
	public myVector U(myVector v, float d, myVector u){myVector r = new myVector(v,d,u); return r._normalize(); }
	public myVector Upt(myPoint v){myVector u = new myVector(v); return u._normalize(); }
	public myVector U(myPoint a, myPoint b){myVector u = new myVector(a,b); return u._normalize(); }
	public myVector U(double x, double y, double z) {myVector u = new myVector(x,y,z); return u._normalize();}
	//float versions of unit vectors
	public myVectorf Uf(myVectorf v){myVectorf u = new myVectorf(v); return u._normalize(); }
	public myVectorf Uf(myVectorf v, float d, myVectorf u){myVectorf r = new myVectorf(v,d,u); return r._normalize(); }
	public myVectorf Uptf(myPointf v){myVectorf u = new myVectorf(v); return u._normalize(); }
	public myVectorf Uf(myPointf a, myPointf b){myVectorf u = new myVectorf(a,b); return u._normalize(); }
	public myVectorf Uf(float x, float y, float z) {myVectorf u = new myVectorf(x,y,z); return u._normalize();}

	//draw a single cell - fclr and sclr
	public void drawFluidCell(myVectorf loc, double sz, double scale){//, int[] fc){
		pushMatrix(); pushStyle();
		fill(0,0,(int) (scale*255),10);//isnt properly transparent - consequence of having z==up, alpha isn't working.
//		fill(fc[0],fc[1],fc[2],10);//isnt properly transparent - consequence of having z==up, alpha isn't working.
//		shininess(0.0f);
		translate(loc.x, loc.y, loc.z);		
//		rotate(HALF_PI,1,0,0);
//		box((float)sz);
		sphere((float) (sz * scale));
		popStyle(); popMatrix();
	}
			
	public myVector normToPlane(myPoint A, myPoint B, myPoint C) {return myVector._cross(new myVector(A,B),new myVector(A,C)); };   // normal to triangle (A,B,C), not normalized (proportional to area)

	public void gl_normal(myVector V) {normal((float)V.x,(float)V.y,(float)V.z);}                                          // changes normal for smooth shading
	public void gl_normal(myVectorf V) {normal(V.x,V.y,V.z);}                                          // changes normal for smooth shading
	public void gl_vertex(myPoint P) {vertex((float)P.x,(float)P.y,(float)P.z);}                                           // vertex for shading or drawing
	public void gl_vertex(myPointf P) {vertex(P.x,P.y,P.z);}                                           // vertex for shading or drawing
	public void showVec( myPoint ctr, double len, myVector v){line(ctr.x,ctr.y,ctr.z,ctr.x+(v.x)*len,ctr.y+(v.y)*len,ctr.z+(v.z)*len);}
	public void showVec( myPointf ctr, double len, myVectorf v){line(ctr.x,ctr.y,ctr.z,ctr.x+(v.x)*len,ctr.y+(v.y)*len,ctr.z+(v.z)*len);}
	public void showVec( myPointf ctr, double len, double _x, double _y, double _z){line(ctr.x,ctr.y,ctr.z,ctr.x+(_x)*len,ctr.y+(_y)*len,ctr.z+(_z)*len);}
	public void show(myPoint P, double r,int fclr, int sclr, boolean flat) {//TODO make flat circles for points if flat
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		if(!flat){
			translate((float)P.x,(float)P.y,(float)P.z); 
			sphereDetail(5);
			sphere((float)r);
		} else {
			translate((float)P.x,(float)P.y,0); 
			this.circle(0,0,(float)r,(float)r);				
		}
		popStyle(); popMatrix();} // render sphere of radius r and center P)
	public void show(myPointf P, double r,int fclr, int sclr, boolean flat) {//TODO make flat circles for points if flat
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		if(!flat){
			translate(P.x,P.y,P.z); 
			sphereDetail(5);
			sphere((float)r);
		} else {
			translate(P.x,P.y,0); 
			this.circle(0,0,(float)r,(float)r);				
		}
		popStyle(); popMatrix();} // render sphere of radius r and center P)
	public void show(myPoint P, double r){show(P,r, gui_Black, gui_Black, false);}
	public void show(myPointf P, double r){show(P,r, gui_Black, gui_Black, false);}
	
	public void show(myPoint P, String s) {text(s, (float)P.x, (float)P.y, (float)P.z); } // prints string s in 3D at P
	public void show(myPoint P, double r, int fclr, int sclr, int tclr, String txt) {
		pushMatrix(); pushStyle(); 
		if((fclr!= -1) && (sclr!= -1)){setColorValFill(fclr); setColorValStroke(sclr);}
		sphereDetail(5);
		translate((float)P.x,(float)P.y,(float)P.z); 
		sphere((float)r); 
		setColorValFill(tclr);setColorValStroke(tclr);
		double d = 1.1 * r;
		show(myPoint.ZEROPT, txt, new myVector(d,d,d));
		popStyle(); popMatrix();} // render sphere of radius r and center P)

	public void show(myPoint P, double r, String s, myVector D){show(P,r, gui_Black, gui_Black, false);pushStyle();setColorValFill(gui_Black);show(P,s,D);popStyle();}
	public void show(myPoint P, double r, String s, myVector D, int clr, boolean flat){show(P,r, clr, clr, flat);pushStyle();setColorValFill(clr);show(P,s,D);popStyle();}
	public void show(myPoint P, String s, myVector D) {text(s, (float)(P.x+D.x), (float)(P.y+D.y), (float)(P.z+D.z));  } // prints string s in 3D at P+D
	public void show(myPoint[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	public void showNoClose(myPoint[] ara) {beginShape(); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape();};                     
	public void show(myPoint[] ara, myVector norm) {beginShape();gl_normal(norm); for(int i=0;i<ara.length;++i){gl_vertex(ara[i]);} endShape(CLOSE);};                     
	public void curveVertex(myPoint P) {curveVertex((float)P.x,(float)P.y);};                                           // curveVertex for shading or drawing
	public void curve(myPoint[] ara) {if(ara.length == 0){return;}beginShape(); curveVertex(ara[0]);for(int i=0;i<ara.length;++i){curveVertex(ara[i]);} curveVertex(ara[ara.length-1]);endShape();};                      // volume of tet 

	public boolean intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C, myPoint X) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		myVector EA=new myVector(E,A), AB=new myVector(A,B), AC=new myVector(A,C); 		double t = (float)(myVector._mixProd(EA,AC,AB) / myVector._mixProd(T,AC,AB));		X.set(myPoint._add(E,t,T));		return true;
	}	
	public myPoint intersectPl(myPoint E, myVector T, myPoint A, myPoint B, myPoint C) { // if ray from E along T intersects triangle (A,B,C), return true and set proposal to the intersection point
		myVector EA=new myVector(E,A), AB=new myVector(A,B), AC=new myVector(A,C); 		
		double t = (float)(myVector._mixProd(EA,AC,AB) / myVector._mixProd(T,AC,AB));		
		return (myPoint._add(E,t,T));		
	}	
	// if ray from E along V intersects sphere at C with radius r, return t when intersection occurs
	public double intersectPt(myPoint E, myVector V, myPoint C, double r) { 
		myVector Vce = V(C,E);
		double CEdCE = Vce._dot(Vce), VdV = V._dot(V), VdVce = V._dot(Vce), b = 2 * VdVce, c = CEdCE - (r*r),
				radical = (b*b) - 4 *(VdV) * c;
		if(radical < 0) return -1;
		double t1 = (b + Math.sqrt(radical))/(2*VdV), t2 = (b - Math.sqrt(radical))/(2*VdV);			
		return ((t1 > 0) && (t2 > 0) ? Math.min(t1, t2) : ((t1 < 0 ) ? ((t2 < 0 ) ? -1 : t2) : t1) );
		
	}	
	
	public void rect(float[] a){rect(a[0],a[1],a[2],a[3]);}				//rectangle from array of floats : x, y, w, h

	
/////////////////////		
///color utils
/////////////////////
	public final int  // set more colors using Menu >  Tools > Color Selector
	  black=0xff000000, 
	  white=0xffFFFFFF,
	  red=0xffFF0000, 
	  green=0xff00FF00, 
	  blue=0xff0000FF, 
	  yellow=0xffFFFF00, 
	  cyan=0xff00FFFF, 
	  magenta=0xffFF00FF,
	  grey=0xff818181, 
	  orange=0xffFFA600, 
	  brown=0xffB46005, 
	  metal=0xffB5CCDE, 
	  dgreen=0xff157901;
	//set color based on passed point r= x, g = z, b=y
	public void fillAndShowLineByRBGPt(myPoint p, float x,  float y, float w, float h){
		fill((int)p.x,(int)p.y,(int)p.z);
		stroke((int)p.x,(int)p.y,(int)p.z);
		rect(x,y,w,h);
		//show(p,r,-1);
	}
	
	public myPoint WrldToScreen(myPoint wPt){return new myPoint(screenX((float)wPt.x,(float)wPt.y,(float)wPt.z),screenY((float)wPt.x,(float)wPt.y,(float)wPt.z),screenZ((float)wPt.x,(float)wPt.y,(float)wPt.z));}

	public int[][] triColors = new int[][] {
		{gui_DarkMagenta,gui_DarkBlue,gui_DarkGreen,gui_DarkCyan}, 
		{gui_LightMagenta,gui_LightBlue,gui_LightGreen,gui_TransCyan}};

	public void setMyColorShiny(myVectorf newColor){  setColorVals((int)newColor.x,(int)newColor.y,(int)newColor.z, 255, 255, 255, 255, 80); }
	public void setMyColorShiny(int[] c){  setColorVals(c[0],c[1],c[2],c[3], 255, 255, 255, 80); }
	public void setMyColorCustomShiny(int[] c, float shinyAmt){ setColorVals(c[0],c[1],c[2],c[3], 255, 255, 255, shinyAmt); }
	
	/*
	*  sets the color to shiny white
	*/
	public void setColorShinyWhite(){		setColorVals(255,255,255, 255, 255, 255, 255, 10.0f);}
	public void setColorFlatWhite(){		setColorVals(255,255,255, 255,55, 55, 55, 1.0f);}
	public void setColorClearWhite(){  		setColorVals(222,222,255,30, 222, 150, 255, 120); }//setColorShinyBrown
	
	public void setColorShinyGreen(){		setColorVals(0,100,40, 255, 55, 255, 55, 80.0f);}
	
	public void setColorShinyBlack(){		setColorVals(0,0,0, 255, 255, 255, 255,300.0f);}
	public void setColorFlatBlack(){		setColorVals(0, 0, 0, 255, 77, 77, 77, 8.0f);}
	public void setColorFlatBlack2(){		setColorVals(0, 0, 0, 255, 25, 25, 25, 3.0f);}
	
	public void setColorSnowGlobeBall(){	setColorSnowGlobeBall(40);}
	public void setColorSnowGlobeBall(int alpha){		setColorVals(0,0,0,alpha, 220,150,255,120.0f);}
	
	public void setColorShinyBrown(){		setColorVals(52,36,3,255, 255,150,20,80.0f);}
	public void setColorFlatBrown(){		setColorVals(78, 42, 5, 255, 77, 77, 77, 3.0f);}
	
	public void setColorShinyMaroon(){	setColorVals(128,0,0,255, 200,155,155, 25.0f); }
	
	public void setColorVals(int[] f, int sp_r, int sp_g, int sp_b, float sh){	setColorVals(f[0],f[1],f[2],f[3] , sp_r, sp_g, sp_b,sh);}
	public void setColorVals(int f_r, int f_g, int f_b, int f_a, 
		  					int sp_r, int sp_g, int sp_b, 
		  					float sh){		  
	    fill(f_r, f_g, f_b, f_a);
	    specular(sp_r, sp_g, sp_b);
	    shininess(sh);
	}//setColorVals
	
	public void setShColorVals(PShape sh, int f_r, int f_g, int f_b, int f_a, 
			int sp_r, int sp_g, int sp_b, 
			float si){		  
		sh.fill(f_r, f_g, f_b, f_a);//(f_r, f_g, f_b, f_a);
		sh.specular(sp_r, sp_g, sp_b);//specular(sp_r, sp_g, sp_b);
		sh.shininess(si);
	}//setColorVals
	public void setClrShinyBlack(PShape sh){		setShColorVals(sh, 0,0,0, 255, 255, 255, 255,300.0f);}
	public void setClrShinyMaroon(PShape sh){		setShColorVals(sh, 128,0,0,255, 200,155,155, 25.0f); }
	public void setClrShinyBrown(PShape sh){		setShColorVals(sh, 52,36,3,255, 255,150,20,80.0f);}
	public void setClrFlatWhite(PShape sh){			setShColorVals(sh,255,255,255, 255,55, 55, 55, 1.0f);}
	public void setClrFlatBlack(PShape sh){		setShColorVals(sh,0, 0, 0, 255, 77, 77, 77, 8.0f);}
	public void setClrFlatBlack2(PShape sh){		setShColorVals(sh,0, 0, 0, 255, 25, 25, 25, 3.0f);}
	public void setClrCustomShiny(PShape sh, int[] c, float shinyAmt){ setShColorVals(sh, c[0],c[1],c[2],c[3], 255, 255, 255, shinyAmt); }  
	
	public static final int[] pipeClr = new int[]{225,200,0,255};
	//ugh
	private void setShClr(PShape sh,int clr){
		switch (clr){
		case 0 : {setClrShinyBlack(sh); return;}
		case 1 : {setClrFlatWhite(sh); return;}
		case 2 : {setClrShinyMaroon(sh); return;}
		case 3 : {setClrShinyBrown(sh); return;}
		case 4 : {setClrFlatBlack2(sh); return;}
		case 5 : {setClrCustomShiny(sh, pipeClr,5.0f);return;}//pipe
		case 6 : {setClrFlatBlack(sh); return;}
		}
	}
	//non-vertex derived pshapes require setfill
	public void setClrSetShinyBlack(PShape sh){		setShSetColorVals(sh, 0,0,0, 255, 255, 255, 255,300.0f);}
	public void setClrSetShinyMaroon(PShape sh){	setShSetColorVals(sh, 128,0,0,255, 200,155,155, 25.0f); }
	public void setClrSetShinyBrown(PShape sh){		setShSetColorVals(sh, 52,36,3,255, 255,150,20,80.0f);}
	public void setClrSetFlatWhite(PShape sh){		setShSetColorVals(sh,255,255,255, 255,55, 55, 55, 1.0f);}
	public void setClrSetFlatBlack(PShape sh){		setShSetColorVals(sh,0, 0, 0, 255, 77, 77, 77, 8.0f);}
	public void setClrSetFlatBlack2(PShape sh){		setShSetColorVals(sh,0, 0, 0, 255, 25, 25, 25, 3.0f);}
	
	public void setClrVals(PShape sh, int[] f, int sp_r, int sp_g, int sp_b, float shn){	setShColorVals(sh,f[0],f[1],f[2],f[3] , sp_r, sp_g, sp_b,shn);}
	public void setClrSetVals(PShape sh, int[] f, int sp_r, int sp_g, int sp_b, float shn){	setShSetColorVals(sh,f[0],f[1],f[2],f[3] , sp_r, sp_g, sp_b,shn);}

	public void setClrSetCustomShiny(PShape sh, int[] c, float shinyAmt){ setShSetColorVals(sh, c[0],c[1],c[2],c[3], 255, 255, 255, shinyAmt); }  
	//set<clr type> functions used with shapes specified in pshape ctor (like SPHERE), whereas regular fill or stroke is only for begin/end shapes
	public void setShSetColorVals(PShape sh, int f_r, int f_g, int f_b, int f_a, 
							int sp_r, int sp_g, int sp_b, 
							float si){		  
		sh.setFill(color(f_r, f_g, f_b, f_a));//(f_r, f_g, f_b, f_a);
		sh.setSpecular(color(sp_r, sp_g, sp_b));//specular(sp_r, sp_g, sp_b);
		sh.setShininess(si);
	}//setColorVals

	
	
	public void setFill(int[] clr){setFill(clr,clr[3]);}
	public void setStroke(int[] clr){setStroke(clr,clr[3]);}		
	public void setFill(int[] clr, int alpha){fill(clr[0],clr[1],clr[2], alpha);}
	public void setStroke(int[] clr, int alpha){stroke(clr[0],clr[1],clr[2], alpha);}
    public void setColorValFill(int colorVal){ setColorValFill(colorVal,255);}
	public void setColorValFill(int colorVal, int alpha){
		switch (colorVal){
			case gui_rnd				: { fill(random(255),random(255),random(255),alpha);break;}
	    	case gui_White  			: { fill(255,255,255,alpha);break; }
	    	case gui_Gray   			: { fill(120,120,120,alpha); break;}
	    	case gui_Yellow 			: { fill(255,255,0,alpha);break; }
	    	case gui_Cyan   			: { fill(0,255,255,alpha);  break; }
	    	case gui_Magenta			: { fill(255,0,255,alpha);break; }
	    	case gui_Red    			: { fill(255,0,0,alpha); break; }
	    	case gui_Blue				: { fill(0,0,255,alpha); break; }
	    	case gui_Green				: { fill(0,255,0,alpha);  break; } 
	    	case gui_DarkGray   		: { fill(80,80,80,alpha); break;}
	    	case gui_DarkRed    		: { fill(120,0,0,alpha);break;}
	    	case gui_DarkBlue   		: { fill(0,0,120,alpha); break;}
	    	case gui_DarkGreen  		: { fill(0,120,0,alpha); break;}
	    	case gui_DarkYellow 		: { fill(120,120,0,alpha); break;}
	    	case gui_DarkMagenta		: { fill(120,0,120,alpha); break;}
	    	case gui_DarkCyan   		: { fill(0,120,120,alpha); break;}	   
	    	case gui_LightGray   		: { fill(200,200,200,alpha); break;}
	    	case gui_LightRed    		: { fill(255,110,110,alpha); break;}
	    	case gui_LightBlue   		: { fill(110,110,255,alpha); break;}
	    	case gui_LightGreen  		: { fill(110,255,110,alpha); break;}
	    	case gui_LightYellow 		: { fill(255,255,110,alpha); break;}
	    	case gui_LightMagenta		: { fill(255,110,255,alpha); break;}
	    	case gui_LightCyan   		: { fill(110,255,255,alpha); break;}    	
	    	case gui_Black			 	: { fill(0,0,0,alpha);break;}//
	    	case gui_TransBlack  	 	: { fill(0x00010100);  break;}//	have to use hex so that alpha val is not lost    	
	    	case gui_FaintGray 		 	: { fill(77,77,77,alpha/3); break;}
	    	case gui_FaintRed 	 	 	: { fill(110,0,0,alpha/2);  break;}
	    	case gui_FaintBlue 	 	 	: { fill(0,0,110,alpha/2);  break;}
	    	case gui_FaintGreen 	 	: { fill(0,110,0,alpha/2);  break;}
	    	case gui_FaintYellow 	 	: { fill(110,110,0,alpha/2); break;}
	    	case gui_FaintCyan  	 	: { fill(0,110,110,alpha/2); break;}
	    	case gui_FaintMagenta  	 	: { fill(110,0,110,alpha/2); break;}
	    	case gui_TransGray 	 	 	: { fill(120,120,120,alpha/8); break;}//
	    	case gui_TransRed 	 	 	: { fill(255,0,0,alpha/2);  break;}
	    	case gui_TransBlue 	 	 	: { fill(0,0,255,alpha/2);  break;}
	    	case gui_TransGreen 	 	: { fill(0,255,0,alpha/2);  break;}
	    	case gui_TransYellow 	 	: { fill(255,255,0,alpha/2);break;}
	    	case gui_TransCyan  	 	: { fill(0,255,255,alpha/2);break;}
	    	case gui_TransMagenta  	 	: { fill(255,0,255,alpha/2);break;}
	    	case gui_OffWhite			: { fill(248,248,255,alpha);break; }
	    	default         			: { fill(255,255,255,alpha);break;}  	    	
		}//switch	
	}//setcolorValFill
	public void setColorValStroke(int colorVal){ setColorValStroke(colorVal, 255);}
	public void setColorValStroke(int colorVal, int alpha){
		switch (colorVal){
	    	case gui_White  	 	    : { stroke(255,255,255,alpha); break; }
 	    	case gui_Gray   	 	    : { stroke(120,120,120,alpha); break;}
	    	case gui_Yellow      	    : { stroke(255,255,0,alpha); break; }
	    	case gui_Cyan   	 	    : { stroke(0,255,255,alpha); break; }
	    	case gui_Magenta	 	    : { stroke(255,0,255,alpha);  break; }
	    	case gui_Red    	 	    : { stroke(255,120,120,alpha); break; }
	    	case gui_Blue		 	    : { stroke(120,120,255,alpha); break; }
	    	case gui_Green		 	    : { stroke(120,255,120,alpha); break; }
	    	case gui_DarkGray    	    : { stroke(80,80,80,alpha); break; }
	    	case gui_DarkRed     	    : { stroke(120,0,0,alpha); break; }
	    	case gui_DarkBlue    	    : { stroke(0,0,120,alpha); break; }
	    	case gui_DarkGreen   	    : { stroke(0,120,0,alpha); break; }
	    	case gui_DarkYellow  	    : { stroke(120,120,0,alpha); break; }
	    	case gui_DarkMagenta 	    : { stroke(120,0,120,alpha); break; }
	    	case gui_DarkCyan    	    : { stroke(0,120,120,alpha); break; }	   
	    	case gui_LightGray   	    : { stroke(200,200,200,alpha); break;}
	    	case gui_LightRed    	    : { stroke(255,110,110,alpha); break;}
	    	case gui_LightBlue   	    : { stroke(110,110,255,alpha); break;}
	    	case gui_LightGreen  	    : { stroke(110,255,110,alpha); break;}
	    	case gui_LightYellow 	    : { stroke(255,255,110,alpha); break;}
	    	case gui_LightMagenta	    : { stroke(255,110,255,alpha); break;}
	    	case gui_LightCyan   		: { stroke(110,255,255,alpha); break;}		   
	    	case gui_Black				: { stroke(0,0,0,alpha); break;}
	    	case gui_TransBlack  		: { stroke(1,1,1,1); break;}	    	
	    	case gui_FaintGray 			: { stroke(120,120,120,250); break;}
	    	case gui_FaintRed 	 		: { stroke(110,0,0,alpha); break;}
	    	case gui_FaintBlue 	 		: { stroke(0,0,110,alpha); break;}
	    	case gui_FaintGreen 		: { stroke(0,110,0,alpha); break;}
	    	case gui_FaintYellow 		: { stroke(110,110,0,alpha); break;}
	    	case gui_FaintCyan  		: { stroke(0,110,110,alpha); break;}
	    	case gui_FaintMagenta  		: { stroke(110,0,110,alpha); break;}
	    	case gui_TransGray 	 		: { stroke(150,150,150,alpha/4); break;}
	    	case gui_TransRed 	 		: { stroke(255,0,0,alpha/2); break;}
	    	case gui_TransBlue 	 		: { stroke(0,0,255,alpha/2); break;}
	    	case gui_TransGreen 		: { stroke(0,255,0,alpha/2); break;}
	    	case gui_TransYellow 		: { stroke(255,255,0,alpha/2); break;}
	    	case gui_TransCyan  		: { stroke(0,255,255,alpha/2); break;}
	    	case gui_TransMagenta  		: { stroke(255,0,255,alpha/2); break;}
	    	case gui_OffWhite			: { stroke(248,248,255,alpha);break; }
	    	default         			: { stroke(55,55,255,alpha); break; }
		}//switch	
	}//setcolorValStroke	
	
    public void setColorValFillAmb(int colorVal){ setColorValFillAmb(colorVal,255);}
	public void setColorValFillAmb(int colorVal, int alpha){
		switch (colorVal){
			case gui_rnd				: { fill(random(255),random(255),random(255),alpha); ambient(120,120,120);break;}
	    	case gui_White  			: { fill(255,255,255,alpha); ambient(255,255,255); break; }
	    	case gui_Gray   			: { fill(120,120,120,alpha); ambient(120,120,120); break;}
	    	case gui_Yellow 			: { fill(255,255,0,alpha); ambient(255,255,0); break; }
	    	case gui_Cyan   			: { fill(0,255,255,alpha); ambient(0,255,alpha); break; }
	    	case gui_Magenta			: { fill(255,0,255,alpha); ambient(255,0,alpha); break; }
	    	case gui_Red    			: { fill(255,0,0,alpha); ambient(255,0,0); break; }
	    	case gui_Blue				: { fill(0,0,255,alpha); ambient(0,0,alpha); break; }
	    	case gui_Green				: { fill(0,255,0,alpha); ambient(0,255,0); break; } 
	    	case gui_DarkGray   		: { fill(80,80,80,alpha); ambient(80,80,80); break;}
	    	case gui_DarkRed    		: { fill(120,0,0,alpha); ambient(120,0,0); break;}
	    	case gui_DarkBlue   		: { fill(0,0,120,alpha); ambient(0,0,120); break;}
	    	case gui_DarkGreen  		: { fill(0,120,0,alpha); ambient(0,120,0); break;}
	    	case gui_DarkYellow 		: { fill(120,120,0,alpha); ambient(120,120,0); break;}
	    	case gui_DarkMagenta		: { fill(120,0,120,alpha); ambient(120,0,120); break;}
	    	case gui_DarkCyan   		: { fill(0,120,120,alpha); ambient(0,120,120); break;}		   
	    	case gui_LightGray   		: { fill(200,200,200,alpha); ambient(200,200,200); break;}
	    	case gui_LightRed    		: { fill(255,110,110,alpha); ambient(255,110,110); break;}
	    	case gui_LightBlue   		: { fill(110,110,255,alpha); ambient(110,110,alpha); break;}
	    	case gui_LightGreen  		: { fill(110,255,110,alpha); ambient(110,255,110); break;}
	    	case gui_LightYellow 		: { fill(255,255,110,alpha); ambient(255,255,110); break;}
	    	case gui_LightMagenta		: { fill(255,110,255,alpha); ambient(255,110,alpha); break;}
	    	case gui_LightCyan   		: { fill(110,255,255,alpha); ambient(110,255,alpha); break;}	    	
	    	case gui_Black			 	: { fill(0,0,0,alpha); ambient(0,0,0); break;}//
	    	case gui_TransBlack  	 	: { fill(0x00010100); ambient(0,0,0); break;}//	have to use hex so that alpha val is not lost    	
	    	case gui_FaintGray 		 	: { fill(77,77,77,alpha/3); ambient(77,77,77); break;}//
	    	case gui_FaintRed 	 	 	: { fill(110,0,0,alpha/2); ambient(110,0,0); break;}//
	    	case gui_FaintBlue 	 	 	: { fill(0,0,110,alpha/2); ambient(0,0,110); break;}//
	    	case gui_FaintGreen 	 	: { fill(0,110,0,alpha/2); ambient(0,110,0); break;}//
	    	case gui_FaintYellow 	 	: { fill(110,110,0,alpha/2); ambient(110,110,0); break;}//
	    	case gui_FaintCyan  	 	: { fill(0,110,110,alpha/2); ambient(0,110,110); break;}//
	    	case gui_FaintMagenta  	 	: { fill(110,0,110,alpha/2); ambient(110,0,110); break;}//
	    	case gui_TransGray 	 	 	: { fill(120,120,120,alpha/8); ambient(120,120,120); break;}//
	    	case gui_TransRed 	 	 	: { fill(255,0,0,alpha/2); ambient(255,0,0); break;}//
	    	case gui_TransBlue 	 	 	: { fill(0,0,255,alpha/2); ambient(0,0,alpha); break;}//
	    	case gui_TransGreen 	 	: { fill(0,255,0,alpha/2); ambient(0,255,0); break;}//
	    	case gui_TransYellow 	 	: { fill(255,255,0,alpha/2); ambient(255,255,0); break;}//
	    	case gui_TransCyan  	 	: { fill(0,255,255,alpha/2); ambient(0,255,alpha); break;}//
	    	case gui_TransMagenta  	 	: { fill(255,0,255,alpha/2); ambient(255,0,alpha); break;}//   	
	    	case gui_OffWhite			: { fill(248,248,255,alpha);ambient(248,248,255); break; }
	    	default         			: { fill(255,255,255,alpha); ambient(255,255,alpha); break; }	    
	    	
		}//switch	
	}//setcolorValFill
	
	//returns one of 30 predefined colors as an array (to support alpha)
	public int[] getClr(int colorVal){		return getClr(colorVal, 255);	}//getClr
	public int[] getClr(int colorVal, int alpha){
		switch (colorVal){
		case gui_rnd					 : { return new int[] {(int) random(255),(int) random(255),(int) random(255),alpha};}
		case gui_LightRnd			     : { return new int[] {(int) random(120)+125,(int) random(120)+125,(int) random(120)+125,alpha};}
		case gui_DarkRnd			     : { return new int[] {(int) random(120)+40,(int) random(120)+40,(int) random(120)+40,alpha};}
    	case gui_Gray   		         : { return new int[] {120,120,120,alpha}; }
    	case gui_White  		         : { return new int[] {255,255,255,alpha}; }
    	case gui_Yellow 		         : { return new int[] {255,255,0,alpha}; }
    	case gui_Cyan   		         : { return new int[] {0,255,255,alpha};} 
    	case gui_Magenta		         : { return new int[] {255,0,255,alpha};}  
    	case gui_Red    		         : { return new int[] {255,0,0,alpha};} 
    	case gui_Blue			         : { return new int[] {0,0,255,alpha};}
    	case gui_Green			         : { return new int[] {0,255,0,alpha};}  
    	case gui_DarkGray   	         : { return new int[] {80,80,80,alpha};}
    	case gui_DarkRed    	         : { return new int[] {120,0,0,alpha};}
    	case gui_DarkBlue  	 	         : { return new int[] {0,0,120,alpha};}
    	case gui_DarkGreen  	         : { return new int[] {0,120,0,alpha};}
    	case gui_DarkYellow 	         : { return new int[] {120,120,0,alpha};}
    	case gui_DarkMagenta	         : { return new int[] {120,0,120,alpha};}
    	case gui_DarkCyan   	         : { return new int[] {0,120,120,alpha};}	   
    	case gui_LightGray   	         : { return new int[] {200,200,200,alpha};}
    	case gui_LightRed    	         : { return new int[] {255,110,110,alpha};}
    	case gui_LightBlue   	         : { return new int[] {110,110,255,alpha};}
    	case gui_LightGreen  	         : { return new int[] {110,255,110,alpha};}
    	case gui_LightYellow 	         : { return new int[] {255,255,110,alpha};}
    	case gui_LightMagenta	         : { return new int[] {255,110,255,alpha};}
    	case gui_LightCyan   	         : { return new int[] {110,255,255,alpha};}
    	case gui_Black			         : { return new int[] {0,0,0,alpha};}
    	case gui_FaintGray 		         : { return new int[] {110,110,110,alpha};}
    	case gui_FaintRed 	 	         : { return new int[] {110,0,0,alpha};}
    	case gui_FaintBlue 	 	         : { return new int[] {0,0,110,alpha};}
    	case gui_FaintGreen 	         : { return new int[] {0,110,0,alpha};}
    	case gui_FaintYellow 	         : { return new int[] {110,110,0,alpha};}
    	case gui_FaintCyan  	         : { return new int[] {0,110,110,alpha};}
    	case gui_FaintMagenta  	         : { return new int[] {110,0,110,alpha};}    	
    	case gui_TransBlack  	         : { return new int[] {1,1,1,alpha/2};}  	
    	case gui_TransGray  	         : { return new int[] {110,110,110,alpha/2};}
    	case gui_TransLtGray  	         : { return new int[] {180,180,180,alpha/2};}
    	case gui_TransRed  	         	 : { return new int[] {110,0,0,alpha/2};}
    	case gui_TransBlue  	         : { return new int[] {0,0,110,alpha/2};}
    	case gui_TransGreen  	         : { return new int[] {0,110,0,alpha/2};}
    	case gui_TransYellow  	         : { return new int[] {110,110,0,alpha/2};}
    	case gui_TransCyan  	         : { return new int[] {0,110,110,alpha/2};}
    	case gui_TransMagenta  	         : { return new int[] {110,0,110,alpha/2};}	
    	case gui_TransWhite  	         : { return new int[] {220,220,220,alpha/2};}	
    	case gui_OffWhite				 : { return new int[] {255,255,235,alpha};}
    	default         		         : { return new int[] {255,255,255,alpha};}    
		}//switch
	}//getClr
	
	public int getRndClrFlagInt(){return (int)random(0,23);}		//return a random color flag value from below
	public int[] getRndClr(int alpha){return new int[]{(int)random(0,255),(int)random(0,255),(int)random(0,255),alpha};	}
	public int[] getRndClr(){return getRndClr(255);	}		
	public Integer[] getClrMorph(int a, int b, double t){return getClrMorph(getClr(a), getClr(b), t);}    
	public Integer[] getClrMorph(int[] a, int[] b, double t){
		if(t==0){return new Integer[]{a[0],a[1],a[2],a[3]};} else if(t==1){return new Integer[]{b[0],b[1],b[2],b[3]};}
		return new Integer[]{(int)(((1.0f-t)*a[0])+t*b[0]),(int)(((1.0f-t)*a[1])+t*b[1]),(int)(((1.0f-t)*a[2])+t*b[2]),(int)(((1.0f-t)*a[3])+t*b[3])};
	}
	
	//used to generate random color
	public static final int gui_rnd = -1;
	//color indexes
	public static final int gui_Black 	= 0;
	public static final int gui_White 	= 1;	
	public static final int gui_Gray 	= 2;
	
	public static final int gui_Red 	= 3;
	public static final int gui_Blue 	= 4;
	public static final int gui_Green 	= 5;
	public static final int gui_Yellow 	= 6;
	public static final int gui_Cyan 	= 7;
	public static final int gui_Magenta = 8;
	
	public static final int gui_LightRed = 9;
	public static final int gui_LightBlue = 10;
	public static final int gui_LightGreen = 11;
	public static final int gui_LightYellow = 12;
	public static final int gui_LightCyan = 13;
	public static final int gui_LightMagenta = 14;
	public static final int gui_LightGray = 15;

	public static final int gui_DarkCyan = 16;
	public static final int gui_DarkYellow = 17;
	public static final int gui_DarkGreen = 18;
	public static final int gui_DarkBlue = 19;
	public static final int gui_DarkRed = 20;
	public static final int gui_DarkGray = 21;
	public static final int gui_DarkMagenta = 22;
	
	public static final int gui_FaintGray = 23;
	public static final int gui_FaintRed = 24;
	public static final int gui_FaintBlue = 25;
	public static final int gui_FaintGreen = 26;
	public static final int gui_FaintYellow = 27;
	public static final int gui_FaintCyan = 28;
	public static final int gui_FaintMagenta = 29;
	
	public static final int gui_TransBlack = 30;
	public static final int gui_TransGray = 31;
	public static final int gui_TransMagenta = 32;	
	public static final int gui_TransLtGray = 33;
	public static final int gui_TransRed = 34;
	public static final int gui_TransBlue = 35;
	public static final int gui_TransGreen = 36;
	public static final int gui_TransYellow = 37;
	public static final int gui_TransCyan = 38;	
	public static final int gui_TransWhite = 39;	
	public static final int gui_OffWhite = 40;
	public static final int gui_LightRnd = 41;
	public static final int gui_DarkRnd = 42;
	
	public static void main(String[] passedArgs) {		
		String[] appletArgs = new String[] { "SnowGlobe2Pkg.SnowGlobeWin" };
		if (passedArgs != null) {   	PApplet.main(PApplet.concat(appletArgs, passedArgs));    } else {   	PApplet.main(appletArgs);    }
	}//main

}
