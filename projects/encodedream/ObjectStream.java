package projects.encodedream;
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
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.concurrent.Semaphore;

public class ObjectStream<E> implements Iterator<E>{
	
	
	private LinkedList<E> queue;
	private boolean closed;
	private Semaphore limit;
	
	public ObjectStream(){
		this(Integer.MAX_VALUE);
	}
	
	public ObjectStream(int limit){
		closed = false;
		queue = new LinkedList<>();
		this.limit = new Semaphore(limit, true);
	}
	
	public ObjectStream(String file, Class<E> clazz) throws IOException, ReflectiveOperationException {
		this(Integer.MAX_VALUE);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String str = null;
		while( (str = reader.readLine()) != null){
			E el = clazz.getConstructor(String.class).newInstance(str);
			add(el);
		}
		reader.close();
		this.close();
	}
	
	
	public void add(E pile){
		if(closed){
			throw new RuntimeException("closed");
		}

		/*while(queue.size()>limit){

			try {
				System.out.println("add wait "+queue.size());
				synchronized(queue){
					queue.notify();
					queue.wait();
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}*/
		try {
			synchronized(queue){
				queue.notify();
			}
//			System.out.println("add wait "+queue.size());
			limit.acquire();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//		System.out.println("add add");
		synchronized(queue){
			queue.addLast(pile);
			queue.notify();
		}
	}
	
	@Override
	public boolean hasNext() {
		if(closed){
			synchronized(queue){
				return queue.size()>0;
			}
		}
		synchronized(queue){
			while(!closed && queue.size() == 0){
				try {
//					System.out.println("hasNext wait");

//					queue.notify();
					queue.wait();


				} catch (InterruptedException e) {
					throw new RuntimeException(e);
				}
			}
		}
//		System.out.println("hasNext return");
		synchronized(queue){
			return queue.size()>0;
		}
	}

	@Override
	public E next() {
		E pile = null;
		limit.release();
		synchronized(queue){
			pile = queue.removeFirst();
			queue.notify();
		}
		return pile;
	}
	
	public void close(){
		closed = true;
		synchronized(queue){
			queue.notify();
		}
	}
	
	
	public void print(PrintStream out){
		while(hasNext()){
			Object o = next();
			out.println(o.toString());
		}
	}
	
}