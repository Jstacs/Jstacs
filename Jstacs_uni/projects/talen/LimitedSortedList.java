/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.talen;

import java.util.Arrays;

import de.jstacs.io.ArrayHandler;
import de.jstacs.utils.ComparableElement;


public class LimitedSortedList<T> {

	private ComparableElement<T, Double>[] list;
	private int curr;
	private boolean expand;
	
	public LimitedSortedList(int limit){
		if(limit > 0){
			list = new ComparableElement[limit];
			expand = false;
		}else{
			list = new ComparableElement[-limit];
			expand = true;
		}
		curr = 0;
	}
	
	public void clear(){
		curr = 0;
	}
	
	public boolean checkInsert(double val){
		return curr < list.length || val > list[0].getWeight();
	}
	
	public ComparableElement<T, Double> getElementAt(int idx){
		if(idx >= curr){
			throw new ArrayIndexOutOfBoundsException();
		}
		return list[idx];
	}
	
	public T getBestElement(){
		if(curr == list.length){
			return list[list.length-1].getElement();
		}else{
			double max = Double.NEGATIVE_INFINITY;
			T maxEl = null;
			for(int i=0;i<curr;i++){
				double temp = list[i].getWeight();
				if(temp > max){
					max = temp;
					maxEl = list[i].getElement();
				}
			}
			return maxEl;
		}
	}
	
	public T getWorstElement(){
		if(curr == list.length){
			return list[0].getElement();
		}else{
			double min = Double.POSITIVE_INFINITY;
			T minEl = null;
			for(int i=0;i<curr;i++){
				double temp = list[i].getWeight();
				if(temp < min){
					min = temp;
					minEl = list[i].getElement();
				}
			}
			return minEl;
		}
	}
	
	public double getBestScore(){
		if(curr == list.length){
			return list[list.length-1].getWeight();
		}else{
			double max = Double.NEGATIVE_INFINITY;
			for(int i=0;i<curr;i++){
				double temp = list[i].getWeight();
				if(temp > max){
					max = temp;
				}
			}
			return max;
		}
	}
	
	public double getWorstScore(){
		if(curr == list.length){
			return list[0].getWeight();
		}else{
			double min = Double.POSITIVE_INFINITY;
			for(int i=0;i<curr;i++){
				double temp = list[i].getWeight();
				if(temp < min){
					min = temp;
				}
			}
			return min;
		}
	}
	
	public void insertAll(LimitedSortedList<T> list2){
		for(int i=0;i<list2.curr;i++){
			insert(list2.list[i].getWeight(),list2.list[i].getElement());
		}
	}
	
	public boolean insert(double val, T element){
		//System.out.println("i: "+val);
		if(curr < list.length){
			list[curr] = new ComparableElement<T, Double>( element, val );
			curr++;
			if(curr == list.length){
				if(expand){
					ComparableElement<T, Double>[] temp = new ComparableElement[(int)(list.length*1.5)];
					System.arraycopy( list, 0,temp, 0, list.length );
					list = temp;
				}else{
					Arrays.sort( list );
				}
			}
			return true;
		}else{
			if(val > list[0].getWeight()){
				ComparableElement<T, Double> el = new ComparableElement<T, Double>( element, val );
				int idx = Arrays.binarySearch( list, el );
				if(idx < 0){
					idx = -idx-1-1;
				}
				System.arraycopy( list, 1, list, 0, idx );
				list[idx] = el;
				//System.out.println(idx);
				
				
				return true;
			}else{
				return false;
			}
		}
	}
	
	public int getLength(){
		return curr;
	}
	
	public ComparableElement<T, Double>[] getSortedList(){
		if(curr == list.length){
			return list.clone();
		}else{
			ComparableElement<T, Double>[] temp = new ComparableElement[curr];
			System.arraycopy( list, 0, temp, 0, curr );
			Arrays.sort( temp );
			return temp;
		}
	}
	
	public String toString(){
		StringBuffer str = new StringBuffer();
		str.append( "[" );
		for(int i=0;i<list.length && (curr == list.length || i < curr);i++){
			str.append( list[i].getWeight()+": "+list[i].getElement().toString()+"; " );
		}
		str.append( "]" );
		return str.toString();
	}
	
}
