package SnowGlobe2Pkg;

import processing.core.PConstants;
import processing.core.PShape;

/**
*  project 2
*  john turner
*  3d scene animation
*  files : 
*  project2Main    (main file, holds setup and draw)
*  project2Global  (holds global contants, methods and variables)
*  snoBall         (defines a snoball object thrown by the snomen at each other
*  snoMan          (defines a snoman object - a collection of spheres and other graphical elements - that is what the main actors are
*  stackPresents   (this file - describes a present stack object - act as collidable terrain, stop snoballs and snowman movement
*  this scene is a few Math.randomly generated snomen having fun in the snow.
*/

public class stackPresents{
	private SnowGlobeWin pa;
	public int ID;                      //this stack of presents' unique id
	public myVectorf location;            //location in world of presents stack
	public int numPresents;             //the number of boxes in this stack
	public myVectorf[] boxDimensions;     //the dimension of each box
	public int[][] boxColor,ribbonColor;     	//the colors of each box
	public float[] boxOrientation;      //the orientation of each box
	public float[] shinyAmt;            //how shiny a present is
	public int wiggle = 5;                    //amount each box might be off center
	public int bxMult = 7;                    //Math.PIxel multiplier for each box

	
	//to pre-render - make once, instead of every display
	public PShape stack;

  
	public stackPresents(SnowGlobeWin _p, int _ID){
		pa = _p;
		ID = _ID; 
		resetLocation();
		numPresents = (int)(Math.random()*8) + 10;			//# presents in this stack
		boxDimensions = new myVectorf[numPresents];
		boxColor = new int[numPresents][];
		ribbonColor = new int[numPresents][];
		boxOrientation = new float[numPresents];
		shinyAmt = new float[numPresents];
		for (int i = 0; i < numPresents; i++){
			float idVal = (i+1) * (pa.groundRadius/300.0f);
			ribbonColor[i] = pa.getClr(pa.gui_DarkRnd);			
			boxDimensions[i] = new myVectorf(idVal * 1.6f * ((Math.random()*.4f) + .8f) ,idVal * 1.4f * ((Math.random()*.4f) + .8f), idVal * 1.4f * ((Math.random()*.4f) + .8f));
			boxOrientation[i] = (float)(Math.random()*pa.TWO_PI);
			boxColor[i] = pa.getClr(pa.gui_rnd);
			shinyAmt[i] = (float)((Math.random()*14)+2 * 30.0f);
		}//init presents' dimensions
		initPresentStack();
	}//constructor
	private void initPresentStack(){
		//draw with "shape(globeBase);"
		stack = pa.createShape(PConstants.GROUP); 
		stack.addChild(buildStack());		
	}//initPieFloor	
	
	private PShape buildStack(){
		float boxHigh = boxDimensions[numPresents-1].z/2.0f;
		PShape shRes = pa.createShape(PConstants.GROUP);
		myVectorf trans = new myVectorf();

		for(int idx = numPresents-1; idx > 0; idx--){
			trans._add(0,0,boxHigh);
			shRes.addChild(buildPresent(idx, trans));
	  		boxHigh = (boxDimensions[idx].z/2.0f) + (boxDimensions[idx-1].z/2.0f); //from center of bottom box to center of next box
	  	}
		trans._add(0,0,boxHigh);
		shRes.addChild(buildPresent(0, trans)); 
		return shRes;
	}
	
	public void resetLocation(){
		location = new myVectorf(0,0,0);
		if (ID != 0){//originally to hide sphere bug in processing - put first stack in center of ball		
		double rad = ((Math.random() * .7f) + .2f) * pa.groundRadius,
				thet = ((Math.random() * PConstants.TWO_PI));      
			location = new myVectorf(rad*Math.cos(thet) , rad * Math.sin(thet),0);	
		}
	}
	  
	//build a box
	private PShape buildPresent(int idx, myVectorf trans){
		PShape shRes = pa.createShape(PConstants.GROUP);
		float[] bDims = new float[]{boxDimensions[idx].x/2.0f,boxDimensions[idx].y/2.0f,boxDimensions[idx].z};
		float divVal = 5.0f/16.0f;
		shRes.translate(trans.x, trans.y,trans.z);
		//box, with space for ribbon
		float[] boxOrient = new float[]{boxOrientation[idx],0,0,1}, zeroAra = new float[]{0,0,0},
				transAra = new float[]{ divVal*boxDimensions[idx].x, divVal * boxDimensions[idx].y,0};
		shRes.addChild(buildBox(idx, boxColor[idx], bDims,  transAra, boxOrient));
		transAra[0] *= -1;
		shRes.addChild(buildBox(idx, boxColor[idx], bDims,  transAra, boxOrient));
		transAra[1] *= -1;
		shRes.addChild(buildBox(idx, boxColor[idx], bDims,  transAra, boxOrient));
		transAra[0] *= -1;
		shRes.addChild(buildBox(idx, boxColor[idx], bDims,  transAra, boxOrient));
		
		//ribbon
		shRes.addChild(buildBox(idx, ribbonColor[idx], new float[]{boxDimensions[idx].x/8,9*boxDimensions[idx].y/8,boxDimensions[idx].z},  zeroAra, boxOrient));
		shRes.addChild(buildBox(idx, ribbonColor[idx], new float[]{9*boxDimensions[idx].x/8,boxDimensions[idx].y/8,boxDimensions[idx].z},  zeroAra, boxOrient));
		return shRes;
	}//buildPresent	
	private PShape buildBox(int idx, int[] bClr, float[] dims, float[] trans, float[] rotVal){
		PShape sh = pa.createShape(pa.BOX, dims);
		pa.setClrSetCustomShiny(sh, bClr, shinyAmt[idx]);
		sh.setStroke(false);
		sh.rotate(rotVal[0],rotVal[1],rotVal[2],rotVal[3]);
		sh.translate(trans[0],trans[1],trans[2]);
		return sh;		
	}//buildBox
	
	public void drawMe(){
		pa.pushMatrix();pa.pushStyle();
		pa.translate(location.x, location.y, location.z);
		pa.shape(stack);
		pa.popStyle();pa.popMatrix();
	}//method draw me
	
	public float getStackWidth(){	return 1.4f*boxDimensions[numPresents-1].x;	}
  
//determines whether the passed location collides with this stack, in the passed direction, with sz being the radius of the cylinder at the passed location
//colLoc is the location of the collision TODO
	public boolean boxCollide(myVectorf loc, float sz, myVectorf[] colLoc){
//	    if (myVector._dist(mover.location, location) < (mover.buttRadL + 1.4*boxDimensions[numPresents-1].x)){
////	      println("collides! sm : " + mover.ID + " loc " + mover.location + " p loc :  " + location + " dist : " + myVector.dist(mover.location, location));     
//	      return true;
//	    }//check if colision  
//	    
		float dist =  myVectorf._dist(loc, location);
		if(dist > sz + getStackWidth() ){return false;}//no col
		//calc col location from collider
		//colLoc[0] = pa.Wf(sz,pa.Uf(location,loc));
		colLoc[0] = pa.Uf(location,loc);
		return true;
		//return (myVectorf._dist(location, location) < (mover.buttRadL + getStackWidth()));//treats bottom box as sphere instead of box for speed
	}
  
}//class presents
