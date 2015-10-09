package de.jstacs.tools;


public class ProgressUpdater {

	private double last;
	private double curr; 
	
	public ProgressUpdater(){
		this.last = -1;
	}
	
	public void setLast(double last){
		this.last = last;
	}
	
	public void setCurrent(double curr){
		this.curr = curr;
	}
	
	public double getPercentage(){
		if(last > 0){
			return curr/last;
		}else{
			return -1;
		}
	}
	
	public void setIndeterminate(){
		this.last = -1;
	}
	
	public boolean isIndeterminate(){
		return last<0;
	}
	
}
