package SnowGlobe2Pkg;

public class myQuaternion {
	myVector v;
	double x,y,z,w, magn, sqMagn;
	public myQuaternion(){v = new myVector();_pset(v,0);	}	
	public myQuaternion(myVector _v, double _w){v = _v;x=v.x;y=v.y;z=v.z; w = _w;	_mag();}	
	public myQuaternion(float _x, float _y, float _z, float _w) {this(new myVector(_x,_y,_z),_w);	}
	public myQuaternion(double _x, double _y, double _z, double _w) {this(new myVector(_x,_y,_z),_w);	}
	public myQuaternion(myQuaternion a) {this(new myVector(a.v),a.w);	}
	
	//call internally whenever values change to keep vector and x,y,z values synched
	protected void _pset(double _x, double _y, double _z, double _w){v.set(_x,_y,_z);x=v.x;y=v.y;z=v.z; w = _w;_mag();}
	protected void _pset(myVector _v, double _w){v.set(_v);x=v.x;y=v.y;z=v.z; w = _w;_mag();}
	//given axis angle representation, convert to quaternion
	public void setFromAxisAngle(double theta, myVector vec){
		double htht = theta/2.0f, sThetH = Math.sin(htht);
		_pset(myVector._mult(vec,sThetH), Math.cos(htht));
	}
	
	public double _mag(){ this.magn = Math.sqrt(this._SqMag()); return magn; }  
	public double _SqMag(){ this.sqMagn =  ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)+ (this.w*this.w)); return this.sqMagn; }  							//squared magnitude
	public double _dot(myQuaternion q){return ((this.x*q.x) + (this.y*q.y) + (this.z*q.z)+ (this.w*q.w));}
	
	public myQuaternion _normalize(){this._mag();if(magn==0){return this;} _pset(x/magn,y/magn,z/magn,w/magn);return this;}
	public static myQuaternion _normalize(myQuaternion _v){_v._mag(); return new myQuaternion(_v.x/_v.magn,_v.y/_v.magn,_v.z/_v.magn,_v.w/_v.magn);}
	//multiply this quaternion by another quaternion -> this * q
	public myQuaternion _qmult(myQuaternion q){
		myVector t1 = myVector._mult(q.v, w), t2 = myVector._mult(v, q.w), t3 = v._cross(q.v); 
		t1._add(t2);t1._add(t3);
		myQuaternion res = new myQuaternion(t1, (this.w * q.w) - v._dot(q.v) );
		return res;
	}
	//give the conjugate of this quaternion
	public myQuaternion _conj(){return new myQuaternion(myVector._mult(v, -1.0),w);}
	
	//rotate the passed vector by this quaternion -> q*q_v*qstar
	public myVector _rot(myVector _v){
		myQuaternion qnorm = myQuaternion._normalize(this),
				q_v = new myQuaternion(_v, 0),
				res1 = qnorm._qmult(q_v),
				conj = qnorm._conj(), 
				res = res1._qmult(conj);
		return res.v;
	}
	
	//rotate toRotVec by passed angle around rVec
	public static myVector _quatRot(double theta, myVector rVec, myVector toRotVec){
		myQuaternion _tmp = new myQuaternion();
		_tmp.setFromAxisAngle(theta,rVec);			//set from passed angle and axis
		return _tmp._rot(toRotVec);	
	}	
	
	//convert this to axis angle - theta,rx,ry,rz
	public double[] toAxisAngle(){
		myQuaternion tmp = new myQuaternion(this);
		if (tmp.w > 1) {tmp._normalize();} 											// if w>1 acos and sqrt will produce errors, this cant happen if quaternion is normalised
		double thet = 2 * Math.acos(tmp.w), s = Math.sqrt(1-tmp.w*tmp.w); 			// assuming quaternion normalised then w is less than 1, so term always positive.
		if (s < 0.0000001) { return new double[]{thet, 1,0,0};} 					// if s close to zero then direction of axis not important
		else {		return new double[]{thet, tmp.x/s, tmp.y/s, tmp.z/s}; }
	}//asAxisAngle
	
	private static myQuaternion _lerp(myQuaternion qa, myQuaternion qb, double t1, double t2){return new myQuaternion(myVector._add(myVector._mult(qa.v,t1), myVector._mult(qb.v,t2)), (qa.w * t1) + (qb.w*t2));}
	
	public static myQuaternion _slerp(myQuaternion qa, myQuaternion qb, double t){
		// Calculate angle between them.
		double cosHalfTheta = qa._dot(qb);
		// if qa=qb or qa=-qb then theta = 0 and we can return qa
		if (cosHalfTheta >= 1.0){	return new myQuaternion(qa);	}
		double halfTheta = Math.acos(cosHalfTheta);
		double s = Math.sqrt(1.0 - cosHalfTheta*cosHalfTheta);
		if (s < 0.0000001){	return _lerp(qa, qb,.5, .5);}
		return _lerp(qa, qb, Math.sin((1 - t) * halfTheta)/s , Math.sin(t * halfTheta) / s);
	}//_slerp
	
	public String toString(){return "vec:"+v.toStrBrf() + "\tw:"+ String.format("%.4f",w);}	
	public String toStringDbg(){return this.toString()+"\tx:"+x+"\ty:"+y+"\tz:"+z;}	

}//myQuaternion

class myQuaternionf {
	myVectorf v;
	public float x,y,z,w,  magn, sqMagn;
	public myQuaternionf(){v = new myVectorf();_pset(v,0);	}	
	public myQuaternionf(myVectorf _x, float _w){v = _x;x=v.x;y=v.y;z=v.z;	w = _w; _mag();	}	
	public myQuaternionf(float _x, float _y, float _z, float _w) {this(new myVectorf(_x,_y,_z),_w);	}
	public myQuaternionf(double _x, double _y, double _z, double _w) {this(new myVectorf(_x,_y,_z),(float)_w);	}
	public myQuaternionf(myQuaternionf a) {this(new myVectorf(a.v),a.w);	}	
	
	//call internally whenever values change to keep vector and x,y,z values synched
	protected void _pset(float _x, float _y, float _z, float _w){v.set(_x,_y,_z);x=v.x;y=v.y;z=v.z; w = _w;_mag();}
	protected void _pset(myVectorf _v, float _w){v.set(_v);x=v.x;y=v.y;z=v.z; w = _w;_mag();}	
	//given axis angle representation, convert to quaternion
	public void setFromAxisAngle(float theta, myVectorf vec){
		float htht = theta/2.0f, sThetH = (float)(Math.sin(htht));
		_pset(myVectorf._mult(vec,sThetH), (float)(Math.cos(htht)));
	}
	//
	public float _mag(){ this.magn = (float)Math.sqrt(this._SqMag()); return magn; }  
	public float _SqMag(){ this.sqMagn =  ((this.x*this.x) + (this.y*this.y) + (this.z*this.z)+ (this.w*this.w)); return this.sqMagn; }  							//squared magnitude
	public float _dot(myQuaternionf q){return ((this.x*q.x) + (this.y*q.y) + (this.z*q.z)+ (this.w*q.w));}

	public myQuaternionf _normalize(){this._mag();if(magn==0){return this;} _pset(x/magn,y/magn,z/magn,w/magn);return this;}
	public static myQuaternionf _normalize(myQuaternionf _v){_v._mag(); return new myQuaternionf(_v.x/_v.magn,_v.y/_v.magn,_v.z/_v.magn,_v.w/_v.magn);}
	//multiply this quaternion by another quaternion -> this * q
	public myQuaternionf _qmult(myQuaternionf q){
		//w1v2 + w2v1 + v1.cross(v2)
		myVectorf t1 = myVectorf._mult(q.v, w), t2 = myVectorf._mult(v, q.w), t3 = v._cross(q.v);
		t1._add(t2);t1._add(t3);
		myQuaternionf res = new myQuaternionf(t1, (this.w * q.w) - v._dot(q.v) );
		return res;
	}
	//give the conjugate of this quaternion
	public myQuaternionf _conj(){return new myQuaternionf(myVectorf._mult(v, -1.0f),w);}
	
	//rotate the passed vector by this quaternion -> q*q_v*qstar
	public myVectorf _rot(myVectorf _v){		
		myQuaternionf qnorm = myQuaternionf._normalize(this),
				q_v = new myQuaternionf(_v, 0),
				res1 = qnorm._qmult(q_v),
				conj = qnorm._conj(), 
				res = res1._qmult(conj);
		return res.v;
	}
	
	//rotate toRotVec by passed angle around rVec
	public static myVectorf _quatRot(float theta, myVectorf rVec, myVectorf toRotVec){
		myQuaternionf _tmp = new myQuaternionf();
		_tmp.setFromAxisAngle(theta,rVec);
		return _tmp._rot(toRotVec);	
	}
	
	//convert this to axis angle - theta,rx,ry,rz
	public float[] toAxisAngle(){
		myQuaternionf tmp = new myQuaternionf(this);
		if (tmp.w > 1) {tmp._normalize();} 														
		float thet = 2 * (float)Math.acos(tmp.w),s = (float)Math.sqrt(1-tmp.w*tmp.w); 			
		if (s < 0.0000001) { return new float[]{thet, 1,0,0};	}			//with s close to 0, thet doesn't matter 
		else {			return new float[]{thet, tmp.x/s, tmp.y/s, tmp.z/s};   }
	}//asAxisAngle
	
	private static myQuaternionf _lerp(myQuaternionf qa, myQuaternionf qb, float t1, float t2){return new myQuaternionf(myVectorf._add(myVectorf._mult(qa.v,t1), myVectorf._mult(qb.v,t2)), (qa.w * t1) + (qb.w*t2));}
	
	public static myQuaternionf _slerp(myQuaternionf qa, myQuaternionf qb, float t){
		// Calculate angle between them.
		float cosHalfTheta = qa._dot(qb);
		// if qa=qb or qa=-qb then theta = 0 and we can return qa
		if (cosHalfTheta >= 1.0){	return new myQuaternionf(qa);	}
		float halfTheta = (float) Math.acos(cosHalfTheta);
		float s = (float) Math.sqrt(1.0 - cosHalfTheta*cosHalfTheta);
		if (s < 0.0000001){	return _lerp(qa, qb,.5f, .5f);}//avoid div by zero - any t will do since denotes axis
		return _lerp(qa, qb, (float) Math.sin((1 - t) * halfTheta)/s , (float) Math.sin(t * halfTheta) / s);
	}//_slerp
	
	public String toString(){return "vec:"+v.toStrBrf() + "\t w:"+ String.format("%.4f",w);}	
	public String toStringDbg(){return this.toString()+"\tx:"+x+"\ty:"+y+"\tz:"+z;}	
	
}//myQuaternionf

