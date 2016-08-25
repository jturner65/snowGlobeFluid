package SnowGlobe2Pkg;

import java.util.TreeMap;

import processing.core.*;

public class snoMan{
  //snowman body is 3 spheres-  head, torso, butt
	protected static SnowGlobeWin pa;
	protected static mySnowGlobeWin win;
	public final int ID;                  //which snoman this is
	 	
	//boolean state flags
	public int[] stFlags;
	public static final int 
		debug = 0,
		turnToTar = 1,                      //should we turn toward our target?
		isWalking = 2,                     //whether the snowman is walking
		isThrowing = 3,                    //snoman is throwing a snoball
		drawSnoBall = 4,
		snowBallInFlight = 5,              //snowball is now in flight;
		hitTarget = 6;
	public static final int numFlags = 7;

	public static final double mveRateMax = 5;
	public static final int happy = 2,
							mad = 1,
							none = 0;
	
	
	//public PShape[] faces;				//array of possible faces, based on emotions, prebuilt
	public static final int numMoods = 3;
	public int feelings,                  //snowman has 3 feelings - none, mad, happy
				targetID,                  //this snoman's target   
				score,                 //how many times this snoman's hit someone
				hits;                  //how many times this snoman's been hit
  
	public float hRadHalf,hRadQtr,
				smH16, smHO7,smH80, smHO4,
				  hrCPO6, hrSPO6,hrSCPO6,
				  hrCPO7, hrSPO7;
	
	
	public myVectorf heading;			//facing direction of this snoman - keep normalized
	
	public float fwdRollAmt,			//amount the body of the snowman "rolls" forward as he is moving 
				smHeight,                //height of snowman
				headRad,                //radius of snowman's head
				hRad60,
				torsoRadS,               //torso radius (in z direction (shorter)
				torsoRadL,               //long torso radius (in xy plane)
				buttRadS,                //butt radius (in z direction (shorter)
				buttRadL,                //long butt radius (in xy plane)
				moveRate,                //how far the snoman travels each step, a multiplier from 2.0 to 20.0
				snoManSqRad,             //distance from center of scene
				hatAngle,                //the angle of the hat on the head
				shinyAmt;                //how shiny this snoman is
  
	public TreeMap<Double,snoMan> neighbors;	//sorted by distance from this snoman list of refs to all snomen
	public Double[] nbrDists;					//dists to all snomen
	
	
	public myVectorf headCenter,            //center of snowman's head - snow man is built along z axis
				  torsoCenter,           //center of snowman's head - snow man is built along z axis
				  buttCenter,           //center of snowman's butt - snow man is built along z axis
				  mySnoballLocation,
				  meToCtr,				//vector from this snowman's center to the center of the snowglobe
				  location;              //center of snowman's contact on ground in world coords
				  
	public int[] myFillColor,           //snoman's current fill color in rgb values
				  mySpecColor,           //snoman's current specular color
				  myHatColor,            //band around snowman's hat
				  
				  mySnoballColor        //color of snoman's snowball
				  ;    //location of current snowball
	
	
	public float animSpd = 4.0f;			//speed of animation
	
	public float[] clickAra;// moveClick, animClick, sbClick;	//various click values to govern animations 
	public static final int
		clMoveIDX = 0,			//animates movement
		clAnimIDX = 1,			//animates arm
		clSBallIDX = 2;			//animates snowball in flight
	public static final int numAnimClicks = 3;
	//flags to check to increment animation, in order of idxs listed for click ara - MUST HAVE same size as # of anim clicks
	public static final int[] animCkFlgs = new int[]{isWalking, isThrowing, snowBallInFlight};
	public static final double[] animPrgChk = new double[]{0.01,0.01, 0.0};
	 
   	public snoMan(SnowGlobeWin _p, mySnowGlobeWin _win, myVectorf _location, float _smHeight, float _theta, int _ID){
  		pa = _p;
  		win = _win;
  		ID = _ID;
    
  		feelings = none;                                //start off with base face
  		shinyAmt = 15;
  		score = 0;                 
  		hits = 0;                  
  		moveRate = (float)(Math.random()*mveRateMax * pa.groundRadius/300.0f)+1;//1.0 -> mveRateMax+1);
  		initFlags();  		
  		setFlag(isWalking, (Math.random() > .2) ? true : false);      //some snowmen might not be walking at first
  		initClickAra();
  		//use to scale rendered object
  		smHeight = _smHeight;
  				
		location = new myVectorf(_location.x, _location.y, _location.z);
		heading = new myVectorf(pa.cos(_theta),pa.sin(_theta),0);//ccw orientation - as _theta increases, snowman turns ccw
		heading._normalize();
		//up is always z= 1
		
		meToCtr = pa.Vf(win.globe.center,location);
  		neighbors  = new TreeMap<Double, snoMan>();  		
    
  		hatAngle = (float)((1 + Math.random()) * pa.PI/3.0);
    
	    myFillColor = pa.getClr(SnowGlobeWin.gui_LightRnd);	    
	    mySpecColor = pa.getClr(SnowGlobeWin.gui_LightRnd);
  	    myHatColor = pa.getClr(SnowGlobeWin.gui_rnd);
	
	    mySnoballLocation = new myVectorf(0,0,0);
	    //initialize snowman snowball color
	    initSnowBall();
	    //snoManSqRad = (((location.x)*(location.x)) + ((location.y)*(location.y)));
	    initSnoMan();
  		//buildFaces();
  	}//snoman constructor (0)  
   	
   	private void initClickAra(){clickAra = new float[numAnimClicks]; for(int i=0;i<numAnimClicks;++i){clickAra[i] = 0.0f;}}
   	
   	private void initSnoMan(){
	    snoManSqRad = (((location.x)*(location.x)) + ((location.y)*(location.y)));// ,.5);
	    //snoManSqRad = (((location.x)*(location.x)) + ((location.z)*(location.z)));
	    //dimensions of snowman's body
	    headRad = smHeight/9.0f;                  //radius of snowman's head
	    hRad60 = headRad/60.0f;
	    hRadHalf = .5f * headRad;
	    hRadQtr = .5f * hRadHalf;
	    smH16 = smHeight/16.0f;
	    smHO7  = smHeight/7.0f;
	    smH80 = smHeight/80.0f;
	    smHO4 = smHeight/4.0f;
	    torsoRadS = 3*smH16;             //torso radius (in z direction (shorter)
	    torsoRadL = smHeight/6.0f;                 //long torso radius (in xy plane)
	    buttRadS = 3*smH16;              //butt radius (in z direction (shorter)
	    buttRadL = smHeight/4.0f;                  //long butt radius (in xy plane)
	    
	    hrCPO7 = headRad * pa.cPiO7;
	    hrSPO7 = headRad * pa.sPiO7;
	    
	    hrCPO6 = headRad * pa.cPiO6;
	    hrSPO6 = headRad * pa.sPiO6;
	    hrSCPO6 = headRad * pa.scPiO6;   		
  		buttCenter = new myVectorf(location.x,location.y, 2.9*smH16);
  		torsoCenter = new myVectorf(location.x, location.y, 8.1*smH16);
  		headCenter = new myVectorf(location.x, location.y, 12.1*smH16);    
  	}
   	
   	
//	//make sure no collisions when building map of snomen by distance
//	private Double chkPutDistInMap(TreeMap<Double, snoMan> map, Double distSq, snoMan snoMan){
//		snoMan chk = map.put(distSq, snoMan);
//		while(chk != null){
//			//replace chk	if not null
//			map.put(distSq, chk);						
//			distSq *= 1.0000000001;//mod distance some tiny amount
//			chk = map.put(distSq, snoMan);				
//		}
//		return distSq;
//	}//chkDistInMap
//	
   	
   	
   	public void buildNeighborList(){
   		neighbors.clear();
   	}
  
  	/*
  	 *  determines a random color for a snowball and a random opponent
  	 */
  	public void initSnowBall(){
  		targetID = ((int)((Math.random()*(win.snoMenCount-2))+1) + ID) % win.snoMenCount;
  		feelings = none;   
  		//initialize snowman snowball color
  		mySnoballColor = pa.getClr(SnowGlobeWin.gui_rnd); 
  	}
   
  	/**
  	 * checks if location's x and y values are greater than the specified radius of the ground
  	 * @param location location vector to check
  	 * @return whether this is past the edge of the world or not
  	 */
  	public boolean chkOutOfBounds(){
  		boolean result = ((meToCtr.sqMagn + (buttRadL*buttRadL)) >= SnowGlobeWin.gSqRad);// (pa.groundRadius * .95)) ;//? true : false;
  		//if (ID == 3){pa.pr(result +"|"+ snoManSqRad +"|"+ buttRadL +"|"+ location.toStrBrf());}
  		if (result){
  			reflectHeading(meToCtr);
  		}
    return result; 
  	}//outofbounds method
  	
  	//reflect heading around t - t not assumed to be normed
  	public void reflectHeading(myVectorf t){
		myVectorf tmp = myVectorf._normalize(t),
		//reflect heading around tmp
		newHeading = myVectorf._add(myVectorf._mult(tmp,-2.0f*heading._dot(tmp)), heading);
		heading.set(newHeading);
		heading._normalize();
  	}
  	
  	
  	public double sqDistFromEye(myPoint lastEyePos){
  		return lastEyePos._SqrDist(location); 
  	}//distFromEye   
//  	/**
//  	 * determines distance from eye, for purposes of determining detail of spheres
//  	 */
//  	public double distFromEye(float[] lastEyePos){
//  		double result;
// 		result = Math.pow(((lastEyePos[0] - location.x) * (lastEyePos[0] - location.x)) + 
//                ((lastEyePos[1] - location.y) * (lastEyePos[1] - location.y)),.5);  
//  		return result; 
//  	}//distFromEye

  	/**
  	 *  updates the various sphere centers
  	 *  with new data if the snowman moves
  	 */
  	public void updateCenters(){    
  		buttCenter.set(location.x,location.y, 2.9*smH16);
  		torsoCenter.set(location.x,location.y, 8.1*smH16);
  		headCenter.set(location.x, location.y, 12.1*smH16);    
  	}

  	/**
  	 *  checks to see if this snoman is currently throwing a snoball
  	 *  @param animClick the click animating the arm throwing
  	 *  @param sbClick the click animating the snoball transit
  	 */  
  	public void checkThrowing(){
  		//whenever throwing, cycle through animation 1 time - piO2's worth of clicks
  		if (getFlag(isThrowing)){
  			if (clickAra[clAnimIDX] > pa.piO2){
  				setFlag(isThrowing,false);
  				setFlag(drawSnoBall,false);
  				clickAra[clAnimIDX] = 0;
  				//win.clickAra[ID+pa.snoMenCount] = 0;
  			}//if click > piO2
  		}//if throwing  
  		else if (!getFlag(snowBallInFlight)){//if not throwing, see if going to throw - don't throw if a snowball in flight already
  			boolean chkVal = (Math.random() < .05);
  			setFlag(isThrowing, chkVal);
  			if (chkVal){
  				setFlag(turnToTar,true);
  				initSnowBall();
  				if (ID == targetID){            System.out.println("this : " + ID + " target : " + targetID);         }//if this id == target id : sanity check
  			}//if throwing - see if we just reset the throwing value
  		}//if throwing/not throwing  
  	}//checkThrowing method
  
  	/**
  	 *  check if this snoman's snoball is in flight
  	 *  if so, check if it hit the target yet, if so, clear out its click and set feelings of thrower and target, as well as target's color
  	 */
  	public void checkSnoballFlying(){  
  		if (getFlag(snowBallInFlight)){
  			if (hitTarget()){
  				//address thrower
  				feelings = happy;
  				score++;
  				setFlag(snowBallInFlight,false);
  				clickAra[clSBallIDX] = 0;
  				//win.clickAra[ID+2*pa.snoMenCount] = 0;
  				//address target
  				win.players[targetID].setHit(mySnoballColor);
  			}
  		}//snoball in flight
  	}//checkSnoballFlying method
  	
  	public void setHit(int[] _ballClr){
		hits++;
		feelings = mad;
		System.arraycopy(_ballClr, 0, myFillColor, 0, _ballClr.length);
  	}
  	
  	
  	//update individual anim clicks based on state of snowman
  	private void updateClicks(){
  		for(int i=0;i<numAnimClicks;++i){
  			incrClick(i, getFlag(animCkFlgs[i]),true);  	  			
  		}
  	}
  	
  	private void incrClick(int idx, boolean chk, boolean doneAnim){if(chk || doneAnim){clickAra[idx]+= mySnowGlobeWin.clickIncr;if (clickAra[idx] > win.clickMax){clickAra[idx] = 0; }	}}

  	/**
  	 *  rotate snowman's orientation
  	 *  @param click used for walking only
  	 *  @param animClick used to animate snowman movements, like throwing
  	 */
  	//moveClick, animClick, sbClick
  	//from double click, double animClick, double sbClick
  	//call rotateAndDrawMe(clickAra[i], clickAra[i+pa.snoMenCount], clickAra[i+(2*pa.snoMenCount)]);
  	public void simulation(){
  		//update vector from center of snowman on ground to center of globe on ground
  		meToCtr.set(win.globe.center,location);// (((location.x) * (location.x)) + ((location.y) * (location.y)));
  		//determines how much the snowman moves forward and how much his body segments "roll" while he's walking
  		float moveRotAmt = (float)(pa.groundRadius/100.0f)*(1-pa.abs(pa.sin(2*clickAra[clMoveIDX] + pa.piO2)));
  		fwdRollAmt = moveRotAmt * buttRadL;
  		//randomly decide whether snoman wants to keep walking or not - only a 5% chance that he'll change his behavior every cycle    
  		toggleFlag(isWalking, .05);
    
//  		//set color to snoman's color
//  		pa.setMyColorCustomShiny(myFillColor,shinyAmt);  
	    //determine if the snoman has finished throwing, and if so turn off throw animation click    
	    checkThrowing();
	    //determine if snoball is in flight, and if so, check if it hit the target.  
	    checkSnoballFlying();

    	//check if snowman is about to fall off world - if so flip heading around normal
    	boolean flippedByFence = chkOutOfBounds();
    	if (!flippedByFence){
	    	if (getFlag(isThrowing) && getFlag(turnToTar)){
	    		setFlag(turnToTar,false);//only turn 1 time per throw
	    		myVectorf newHead = new myVectorf(location,win.players[targetID].location);
	    		newHead._normalize();
	    		heading = newHead;
	    	}//if throwing and need to face target
	    	else if (getFlag(isWalking)){ //if not flipped by fence and walking, then check 
	    		if ((Math.random() > .750) && (moveRotAmt < 0.01)){      heading = myQuaternionf._quatRot((float)(Math.random()*pa.TWO_PI) - pa.PI, new myVectorf(0,0,-1), heading); heading._normalize();   }
		    	//check if snowman is about to walk into another snowman - if so, turn him around by 180 degrees
		    	if (win.collidingSM(this) != -1){//also flips heading
		    		setFlag(isWalking, (Math.random()*100) > 5) ;       
		    	}
	    	} 
    	}
    	if (getFlag(isWalking) && (moveRotAmt >0.01)){ 		  moveMe(moveRotAmt);    }//if walking
	    	//rotate this snowman to the correct orientation
    
	 //   drawMe(moveRotAmt);
  	}//rotateAndDrawMe
  
  	/**
  	 *  returns whether the snoball hit the target snoman yet
  	 *  @param whether a hit is registered
  	 */
  	public boolean hitTarget(){
  		return (( Math.abs(win.snoballs[ID].location.x - win.players[targetID].location.x) 
  	             + Math.abs(win.snoballs[ID].location.y - win.players[targetID].location.y)) < (win.players[targetID].buttRadL));   
  	}
  
  
  	/**
  	 *  modifies the snoman's location on the ground based on any movement, if happening, and then
  	 *  updates each component's set of axes,  making sure not to have him move backward in the cyclic motion of movement
  	 */
  	public void moveMe(float movRotAmt){
  		location.x = location.x + (moveRate*movRotAmt * heading.x);
  		location.y = location.y - (moveRate*movRotAmt * heading.y);
  		updateCenters(); 
  	}//moveMe
  	
	private void initFlags(){stFlags = new int[1 + numFlags/32]; for(int i = 0; i<numFlags; ++i){setFlag(i,false);}}
	public void setFlag(int idx, boolean val){
		int flIDX = idx/32, mask = 1<<(idx%32);
		stFlags[flIDX] = (val ?  stFlags[flIDX] | mask : stFlags[flIDX] & ~mask);
		switch (idx) {//special actions for each flag
			case debug : {break;}
			case turnToTar: {break;}
			case isWalking : {break;}
			case isThrowing : {break;}
			case drawSnoBall: {break;}
			case snowBallInFlight: {break;}
			case hitTarget: {break;}
		}
	}//setFlag	
	public boolean getFlag(int idx){int bitLoc = 1<<(idx%32);return (stFlags[idx/32] & bitLoc) == bitLoc;}	
	//toggle value of passed flag with given probability
	public void toggleFlag(int idx, double prob){
		//walking = (Math.random() < .950) ? walking : !walking;
		if(Math.random() <= prob){setFlag(idx,!getFlag(idx));}		
	}
  
  	/**
  	 *  draw's the current snowman at his current center
  	 *  @param lastEyePos double ara of x,y,z coords of eye position
  	 *  @param click timer for this particular snowMan
  	 */
  	public void drawMe(){    
	    pa.pushMatrix();pa.pushStyle();
	    //first move axes to center of snowman - all subsequent translations keep this value
  		//set color to snoman's color
  		pa.setMyColorCustomShiny(myFillColor,shinyAmt);  
    	pa.translate(location.x, location.y, location.z);
		//determine how detailed to draw the speheres of this snowman based on his distance from the camera
  		//double dist2Eye = sqDistFromEye(pa.c.getEyeLoc());
  		int sphereDet = 16;// +  Math.max(0,(int)(10-(dist2Eye/20)));
  		pa.sphereDetail(sphereDet, sphereDet+1);
  		float heading_thet = myVectorf._angleBetween(heading, new myVectorf(1,0,0));
  		pa.rotate(heading_thet,0,0,1);
  		//draw snowman's butt
  		pa.pushMatrix();
  			pa.translate(0,0,buttCenter.z);
	    	pa.scale((buttRadL/buttRadS), (buttRadL/buttRadS), 1);      //need to "squish" butt into elongated shape
	    	pa.rotate(((pa.piO6) * fwdRollAmt),0,0,-1);
	    	pa.sphere(buttRadS);       
	    pa.popMatrix();
    
	    //the rest of the snowman will move together based on whether he's walking or not
	    //should move forward by a certain amount for a "step" and then come back to center 
	    //
	    //draw snowman's torso
	    pa.pushMatrix();
	    	pa.translate(fwdRollAmt, 0, (torsoCenter.z));
	    	pa.scale((torsoRadL/torsoRadS), (torsoRadL/torsoRadS), 1);
	    	pa.sphere((torsoRadS)); 
	    pa.popMatrix();
        
      //draw snowman's head 
      	pa.pushMatrix();
      		pa.translate(fwdRollAmt,0,(headCenter.z));
      		pa.sphere((headRad)); 
      	pa.popMatrix();//head
      	
	    //draw snowman's right arm- may be throwing a snowball
	    pa.pushMatrix();
	    	pa.translate(fwdRollAmt, (torsoRadL), (torsoCenter.z));
	    	pa.rotate((pa.PI/3.0f),-1,0,0);
	    	if (getFlag(isThrowing)) {
	    		drawThrow(smH80,smHO4, 4);
	    		if (getFlag(drawSnoBall)){//draw snoball
	    			mySnoballLocation.set(location.x+(1.5f*buttRadL*heading.y),location.y+(1.5f*buttRadL*heading.x), headCenter.z);// = new myVector(modelX(0,0,0), modelY(0,0,0), modelZ(0,0,0));
	    			win.snoballs[ID].setColor(mySnoballColor);
	    			win.snoballs[ID].location.set(mySnoballLocation);
	    			win.snoballs[ID].drawMe(clickAra[clSBallIDX]);
	    		}//if draw snoball
	    	} else {		drawArm(smH80,smHO4, 4, false); 	}
        pa.popMatrix();    
      
        //draw snowman's left arm
        pa.pushMatrix();
        	pa.translate(fwdRollAmt, (-torsoRadL),(torsoCenter.z));
      		pa.rotate((-pa.PI/3.0f),-1,0,0);
      		drawArm(smH80, -smHO4, 4,true);  
      	pa.popMatrix();

      	//draw snowman's face
      	pa.pushMatrix();
      		pa.translate(fwdRollAmt,0, headCenter.z);
      		drawFace();
      	pa.popMatrix();
       	pa.popStyle();pa.popMatrix();//pop matrix for translation and rotation to heading

       	
    	//animate flying snowball
	    if (getFlag(snoMan.snowBallInFlight)) {   	win.snoballs[ID].drawMe(clickAra[clSBallIDX]);    }
    	updateClicks();
  	}//draw me
  	
    
  	/**
  	 *  draw proper face based on this snoman's emotional state of mind 
  	 */
  	public void drawFace(){
  		drawEyes();
  		drawEyebrows(feelings * pa.PI/4.0f); //-1 is happy, 1 is mad, 
    	drawMouth();
 		drawNose();  
 	   	drawTopHat();   
  	}//drawface
  	
  	
 	
  	/**
  	 *  2 eyes made out of coal
  	 *  anthracite, to be specific.  nice, shiny anthracite
  	 */
  	public void drawEyes(){
  		pa.setColorShinyBlack();
  		//draw eyes
  		pa.pushMatrix();    
  			pa.translate(hrCPO7, hrSPO7*pa.cPiO6, hrSPO7*pa.cPiO6);
  			pa.sphere(smH80);
  		pa.popMatrix();
  		pa.pushMatrix();
  			pa.translate(hrCPO7, -hrSPO7*pa.cPiO6, hrSPO7*pa.cPiO6);
   			pa.sphere(smH80);
  		pa.popMatrix();//eyebrows
  	}
    
  	/**
  	 *  expressive eyebrows
  	 *  feelings : 0 - none, 1 - mad, -1 - happy
  	 *  rotAmt is how much to rotate based on feelings
  	 */
  	public void drawEyebrows(float rotAmt){
  		//set color to flat black
  		pa.setColorFlatBlack2();
  		//draw eye brows
  		pa.pushMatrix();
  			pa.translate((1.1f*hrCPO6), 
  						(headRad * pa.scPiO6) - hRadQtr, 
  						(headRad * pa.scPiO6 + headRad/10.0));
			pa.translate(0, hRadQtr, 0);
			pa.rotate(rotAmt, 1,0,0);
			pa.translate(0, -hRadQtr, 0);
  			//pa.translate(0,(headRad/20.0),0);
  			pa.drawCylinder(hRad60,hRad60,hRadHalf,4,false,false);
      
  		pa.popMatrix();
  		pa.pushMatrix();  			
  			pa.translate((1.1f*hrCPO6), 
  					(-1*headRad * pa.scPiO6) - hRadQtr, 
  					(headRad * pa.scPiO6 + headRad/10.0));
			pa.translate(0, hRadQtr, 0);
			pa.rotate(-rotAmt, 1,0,0);
			pa.translate(0, -hRadQtr, 0);
  			//pa.translate(0,(-headRad/20.0),0);
  			pa.drawCylinder(hRad60,hRad60,hRadHalf,4,false,false);      
  		pa.popMatrix();//eyebrows    
    
  	}//drawEyebrows
    
  	/**
  	 *  an expressive mouth, holding a pa.PIpe
  	 *  feelings : 0 - none, 1 - mad, -1 - happy
  	 */
 	public void drawMouth(){
  		//draw mouth, based on feelings
  		pa.setColorShinyMaroon();
  		pa.pushMatrix();
			pa.translate(hrCPO6,
	        		(headRad* pa.sPiO12),
	        		(-headRad * pa.sPiO6* pa.cPiO12));
	  		//if (feelings == happy){     
	  		pa.sphere(smH80);         
//	  		}//feelings = happy
//	  		else if (feelings == mad ){     
//	        	pa.sphere(smH80);        
//	  		}//feelings = mad
//	  		else{//if none        
//	  			pa.sphere(smH80);        
//	      	}//feelings = none    
  		pa.popMatrix();
  		drawPipe();        
  	}//drawMouth
  
  	/*
  	 *  carrot nose
  	 */
  	public void drawNose(){
  		//set color to orange for nose
  		pa.setMyColorCustomShiny(new int[]{225,120,0,255}, 4.0f);
  		pa.pushMatrix();
  			pa.translate(headRad,0,0);
  			pa.rotate(-pa.piO2,0,1, 0);
  			//drawAxes();
  			pa.scale(0.1f ,0.1f,0.6f);
  			pa.sphereDetail(4);
  			pa.sphere(headRad);
  		pa.popMatrix();  
  	}//drawNose
  
  	/**
  	 *  corncob pipe
  	 */
  	public void drawPipe(){
    //set color to yellow for pipe
  		pa.pushStyle();
		pa.setMyColorCustomShiny(pa.pipeClr, 5.0f);

		pa.pushMatrix();
			pa.translate(hrCPO6,
						(headRad* pa.sPiO12),
						(-hrSPO6* pa.cPiO12));
			pa.rotate((3*pa.piO2),0,0,1);
			pa.rotate((-pa.piO6),1,0,0);
			pa.drawCylinder(headRad/90.0f,headRad/90.0f,headRad/2.0f,3,false,false);  //pipeestem
			pa.translate(0,(headRad/2.0),0);
			pa.rotate((pa.piO2),1,0,0);
			pa.drawCylinder(headRad/10.0f,headRad/10.0f,headRad/5.0f,7,true,true);  //pipehead
			pa.translate(0,(headRad/4.95f),0);
			pa.rotate((pa.piO2),1,0,0);			
			pa.setColorFlatBlack2();  
			pa.ellipse(0,0,(headRad/5.0f),(headRad/5.0f));                         //black circle for pipe
		pa.popMatrix();
		pa.popStyle();
  	}//drawpa.PIpe

  	/**
  	 *  what's a snoman without his tophat?
  	 */  
  	public void drawTopHat(){//(bottom,top,height, drawbottom, drawtop)
  		pa.pushStyle();
  		pa.setColorFlatBlack();
  		pa.pushMatrix();
  			pa.rotate((hatAngle),1,0,0);
  			pa.translate(0,(headRad * .8),0);
  			//main hat part
  			pa.drawCylinder(3*headRad/5.0f,4*headRad/5.0f, 1.2f*headRad, 16, false,true);  
  			//brim
  			pa.drawCylinder(headRad*1.0f,headRad*1.2f,headRad*.04f, 16,true, true);
  			//band around hat - random color
  			pa.translate(0,(0.1*headRad),0);
  
  			setColorTopHatBand();
  			//band
  			pa.drawCylinder(3.5f*headRad/5.0f, 3.3f*headRad/5.0f,headRad/4.0f, 16, false, false);
  		pa.popMatrix();
  		pa.popStyle();
  	}//drawtophat for this snowman
  
  	public void drawArm(float rad, float armLength, int sides,boolean leftSide) {
  		//set color of arms to be flat brown 
  		pa.setColorFlatBrown();
  		pa.pushStyle();
  		//upper and lower arms are cylinders
  		pa.translate(0,(-armLength/8.0),0);
  		pa.drawCylinder(rad,rad, 5*armLength/8.0f, sides, false, false);  
  		pa.translate(0,(5*armLength/8.0f),0);
  		pa.drawCylinder(rad,rad, armLength/2.0f, sides, false, false);  
    
  		//now need to draw hand - use color of tophat band
  		setColorTopHatBand();
  		pa.pushMatrix();
  			pa.translate(0,(armLength/2.0),0);
  			//drawAxes();
  			if (leftSide){      pa.rotate((pa.PI),1,0,0);    }
  			pa.scale(2,3,1);
  			pa.sphere((headRad/7.0f));
  		pa.popMatrix();    
  		//setColorDefault();
  		pa.popStyle();
     
  	}//draw arm
  
  	/*
  	 * animates throwing a snowball - used instead of drawArm for right arm if snowball is being thrown
   	 */
  	public void drawThrow(float rad, float armLength, int sides){
  		//this is just rotations at the joints - driven by click
  		float clickVal = PApplet.sin(clickAra[clAnimIDX] * animSpd), absClkVal = PApplet.abs(clickVal);

  		//set color of arms to be flat brown 
  		pa.setColorFlatBrown();
  		//upper and lower arms are cylinders
  		//drawAxes();
  		pa.translate(0,(-2.0f*armLength/8.0f),0);
  		//shoulder rotation
  		pa.rotate((2*PApplet.abs(clickVal * .5f)),1,0,0);
    	pa.rotate((-1.2f*absClkVal),0,0,1);
    	pa.translate(0,armLength/8.0f,0);
    	pa.drawCylinder(rad,rad, 5*armLength/8.0f, sides, false, false);  
    	pa.translate(0,(5*armLength/8.0f),0);
    	if (clickVal > .9){    setFlag(drawSnoBall,true); }
    	else if (clickVal < 0) {     
    		clickAra[clAnimIDX] = 0;
    		if (getFlag(drawSnoBall)){
    			setFlag(drawSnoBall,false);
    			setFlag(snowBallInFlight,true);
    		}//launch snoball subroutine
    	} 
    	pa.rotate(-2.5f*absClkVal,0,0,1);
    	pa.drawCylinder(rad,rad, armLength/2.0f, sides, false, false);  
    
    	//now need to draw hand - use color of tophat band
    	setColorTopHatBand();
    	pa.translate(0,armLength/2.0f,0);
    	pa.pushMatrix();
    		pa.scale(2,3,1);
    		pa.sphere((headRad/7.0f));        //draw mitten
    		//pa.translate(1,0,0);
    		//need to draw thumb here   
    	pa.popMatrix();
  	}//drawThrow method  
  
	private void setColorTopHatBand(){pa.setColorVals(myHatColor, 200, 200, 200, 30.0f);}	
	private void setShColorVals(PShape sh, int f_r, int f_g, int f_b, int f_a, 
		  					int sp_r, int sp_g, int sp_b, 
		  					float si){		  
	    sh.fill(f_r, f_g, f_b, f_a);
	    sh.specular(sp_r, sp_g, sp_b);
	    sh.shininess(si);
	}//setColorVals

		
	/**
	* returns relevant info about the snowman
	* @return info about this snowman
	*/
	public String toString(){
	  String info = "height : " + smHeight + " width : " + smHeight/2.0f + " location in world : " + location;
	 
	  return info;
	}

}//snoman

//class holds pre-rendered version of snowman held as pshapes, to speed up drawing
class SMRndrdObj{
	protected static SnowGlobeWin pa;						//papplet
	protected static mySnowGlobeWin win;					//display window
	
	public static final float smHeight = 100.0f;
	
	private float hRadHalf,hRadQtr, smH16, smHO7,smH80, smHO4,
	  hrCPO6, hrSPO6,hrSCPO6, hrCPO7, hrSPO7,
		headRad,                //radius of snowman's head
		hRad60,
		torsoRadS,               //torso radius (in z direction (shorter)
		torsoRadL,               //long torso radius (in xy plane)
		buttRadS,                //butt radius (in z direction (shorter)
		buttRadL;                //long butt radius (in xy plane)

	private myVectorf headCenter,            //center of snowman's head - snow man is built along z axis
	  torsoCenter,           //center of snowman's head - snow man is built along z axis
	  buttCenter;           //center of snowman's butt - snow man is built along z axis

	
	
	public SMRndrdObj(SnowGlobeWin _p, mySnowGlobeWin _win){
		pa=_p;win=_win;
		
		
	}
	
   	private void initSnoMan(){
	    //build everything off of smHeight
	    //snoManSqRad = (((location.x)*(location.x)) + ((location.z)*(location.z)));
	    //dimensions of snowman's body
	    headRad = smHeight/9.0f;                  //radius of snowman's head
	    hRad60 = headRad/60.0f;
	    hRadHalf = .5f * headRad;
	    hRadQtr = .5f * hRadHalf;
	    smH16 = smHeight/16.0f;
	    smHO7  = smHeight/7.0f;
	    smH80 = smHeight/80.0f;
	    smHO4 = smHeight/4.0f;
	    torsoRadS = 3*smH16;             //torso radius (in z direction (shorter)
	    torsoRadL = smHeight/6.0f;                 //long torso radius (in xy plane)
	    buttRadS = 3*smH16;              //butt radius (in z direction (shorter)
	    buttRadL = smHeight/4.0f;                  //long butt radius (in xy plane)
	    
	    hrCPO7 = headRad * pa.cPiO7;
	    hrSPO7 = headRad * pa.sPiO7;
	    
	    hrCPO6 = headRad * pa.cPiO6;
	    hrSPO6 = headRad * pa.sPiO6;
	    hrSCPO6 = headRad * pa.scPiO6;   		
  		buttCenter = new myVectorf(0,0, 2.9*smH16);
  		torsoCenter = new myVectorf(0,0, 8.1*smH16);
  		headCenter = new myVectorf(0,0,12.1*smH16);    
  	}
   	
	
	
	

}
