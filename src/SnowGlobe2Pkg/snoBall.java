package SnowGlobe2Pkg;

import processing.core.PApplet;
import processing.core.PConstants;

public class snoBall{
	//a snoball is a sphere thrown by a snoman at another snoman
	//it has a color which changes the target snoman's color when it hits
	protected static SnowGlobeWin pa;
	protected static mySnowGlobeWin win;
	public final float clickMult = .01f;
	public double ballRad;
	public final int ID;
	public int targetID;
	public snoMan thrower;
	public int[] ballColor;
	public myVectorf location, targetLoc,  start;
	protected final float zTrajMult;
  //PMatrix3D throwMatrix;
  
  	public snoBall(SnowGlobeWin _p,mySnowGlobeWin _win, int _ID){
  		pa = _p; win=_win;
 		ID = _ID;
 		thrower = win.players[ID];
 		ballRad = (thrower.headRad/2.5); 		
   		targetID = ID;
  		location = new myVectorf(0,0,-100);
  		start = new myVectorf(0,0,-100);
  		ballColor = new int[]{150,150,150,255}; 
  		zTrajMult = pa.groundRadius/15.0f;
  	}//snoball constructor
  	
  /**
  *  draws the snoball
  */
  	public void drawMe(float _click){
  		targetID = thrower.targetID;
  		ballColor = thrower.mySnoballColor;
  		//set the ball's color
  		pa.setMyColorShiny(ballColor);
  		pa.sphereDetail(10); 
  		fly(_click);
  		pa.pushMatrix();
  		//don't translate if being translated already in snoman's routine
  		if (!thrower.getFlag(snoMan.drawSnoBall)){pa.translate((float)location.x , (float)location.y , (float)location.z );  }
//      if (ID == 0){
//        println("---in snoball "+ ID +"'s draw : " + location + "target : " + p.players[targetID].location + " click :" + click);
//      }
  			pa.sphere((float)ballRad);
  		pa.popMatrix();    
  	}//snoball drawme method
  	
  	public void setColor(int[] clr){	System.arraycopy(clr, 0, ballColor, 0, clr.length);  	}
  	
  	/**
  	 *  this will compute the current location for the snowball, based on the click and the target's movement.
  	 * 	TODO replace with prediction of trajectory
  	 */
  	public void fly(float _click){
  		float clikMlt = _click*clickMult;
  		location.set(location.x + clikMlt*(win.players[targetID].location.x - location.x) , 
                      location.y + clikMlt*(win.players[targetID].location.y - location.y) , 
                      (thrower.torsoCenter.z + clikMlt*(win.players[targetID].torsoCenter.z - thrower.torsoCenter.z)) + (zTrajMult*PApplet.sin((_click * PConstants.PI)/1.5f)));
   	}
  	public String toString(){
  		String result = "";
  		result += location.x + " | " + location.y + " | " + location.z + " | source id " + ID + " | target id : " + targetID;    
  		return result;
  	}
}//class snoball
