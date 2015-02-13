package de.jstacs.clustering.distances;

import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import projects.motifComp.DeBruijnMotifComparison;
import de.jstacs.data.DeBruijnSequenceGenerator;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.CyclicSequenceAdaptor;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.sequenceScores.statisticalModels.StatisticalModel;


public class SequenceScoreDistance extends DistanceMetric<StatisticalModel> {

	private CyclicSequenceAdaptor[] seqs;
	
	public SequenceScoreDistance(DiscreteAlphabet alphabet, int n) throws OperationNotSupportedException, WrongAlphabetException, WrongSequenceTypeException{
		this.seqs = DeBruijnSequenceGenerator.generate( alphabet, n );
	}
	
	public double[][] getProfile(StatisticalModel o, boolean rc) throws Exception{
		return DeBruijnMotifComparison.getProfilesForMotif( seqs, o, rc );
	}
	
	public double getDistance(double[][] profiles1, double[][] profiles1Rc, double[][] profiles2, int maxShift) throws Exception{
		double fwd = DeBruijnMotifComparison.compare( profiles1, profiles2, maxShift ).getSecondElement();
		double rc = DeBruijnMotifComparison.compare( profiles1Rc, profiles2, maxShift ).getSecondElement();
		
		return 1.0-Math.max( fwd, rc );
	}
	
	public double getDistance(double[][] profiles1, double[][] profiles1Rc, StatisticalModel o2, int motif1Length) throws Exception{
		int maxShift = Math.max( motif1Length - (int)Math.floor(o2.getLength()/3), o2.getLength() - (int)Math.floor( motif1Length/3 ) );
		double[][] profiles2 = getProfile( o2, false );
		return getDistance( profiles1, profiles1Rc, profiles2, maxShift );
	}
	
	@Override
	public double getDistance( StatisticalModel o1, StatisticalModel o2 ) throws Exception {
		
		int motif1Length = o1.getLength();
		
		double[][] profiles1 = getProfile( o1, false );
		double[][] profiles1Rc = getProfile( o1, true );
		
		return getDistance( profiles1, profiles1Rc, o2, motif1Length );
		
	}
	
	public double[][] getPairwiseDistanceMatrix(int numThreads, StatisticalModel... objects) throws Exception {
		if(numThreads == 1){
			double[][] matrix = new double[objects.length][];
			for(int i=0;i<matrix.length;i++){
				matrix[i] = new double[i];
				double[][] prof1 = this.getProfile( objects[i], false );
				double[][] prof1Rc = this.getProfile( objects[i], true );
				fillDistanceRow( i, matrix[i], prof1, prof1Rc, this, objects );
			}
			return matrix;
		}else{
			LinkedList<Worker> available = new LinkedList<Worker>();
			Worker[] workers = new Worker[numThreads-1];
			for(int i=0;i<workers.length;i++){
				workers[i] = new Worker( objects, available, i );
				available.add( workers[i] );
				(new Thread(workers[i])).start();
			}
			
			double[][] matrix = new double[objects.length][];
			for(int i=0;i<matrix.length;i++){
				System.out.println("main in row "+i);
				matrix[i] = new double[i];
				double[][] prof1 = null;
				double[][] prof1Rc = null;
				int length = 0;
				synchronized(objects[i]){
					prof1 = this.getProfile( objects[i], false );
					length = objects[i].getLength();
					prof1Rc = this.getProfile( objects[i], true );
				}
				int a = 0;
				while(a == 0){
					synchronized( available ) {
						a = available.size();
					}
					if(a == 0){
						synchronized(this){
					//		System.out.println("main waits");
							wait(1000);
						}
					}
				}
				Worker curr = null;
				synchronized( available ) {
					curr = available.pop();
				}
			//	System.out.println("main found worker "+curr.index);
				curr.setRow( matrix[i], prof1, prof1Rc, length );
			//	System.out.println("main set row in worker "+curr.index);
			}
			int a = 0;
			while(a < workers.length){
				synchronized( available ) {
					a = available.size();
				}
				if(a == workers.length){
					synchronized(this){
		//				System.out.println("main waits");
						wait(100);
					}
				}
			}
			
			for(int i=0;i<workers.length;i++){
				workers[i].stop();
			}
			return matrix;
		}
	}
	
	private static void fillDistanceRow(int i, double[] row, double[][] prof1, double[][] prof1Rc, SequenceScoreDistance metric, StatisticalModel... objects) throws Exception{
		for(int j=0;j<i;j++){
			row[j] = metric.getDistance( prof1, prof1Rc, objects[j], objects[i].getLength() );
		}
	}

	
	private class Worker implements Runnable{

		private int index;
		private boolean stop;
		private double[] row;
		private double[][] prof1;
		private double[][] prof1Rc;
		private SequenceScoreDistance metric;
		private StatisticalModel[] objects;
		private LinkedList<Worker> available;
		private int length;
		
		private Worker(StatisticalModel[] objects, LinkedList<Worker> available, int index){
			this.objects = objects;
			this.metric = SequenceScoreDistance.this;
			this.available = available;
			this.index = index;
		}
		
		private void setRow(double[] row, double[][] prof1, double[][] prof1Rc, int length){
			//System.out.println("set row");
			synchronized( this ) {
				this.row = row;
				this.prof1 = prof1;
				this.prof1Rc = prof1Rc;
				this.length = length;
				
				//System.out.println("notified this in set");
				notify();
			}
		}
		
		private void stop(){
			this.stop = true;
			synchronized(this){
				//System.out.println("notified this in stop");
				notify();
			}
		}
		
		@Override
		public void run() {
			synchronized(this){
				try {
					//System.out.println("worker waits");
					wait();
				} catch ( InterruptedException e ) {
					e.printStackTrace();
				}
			}
			while(!stop){
				
				if(row != null){
					//System.out.println("worker "+index+" works");
					try {
						fillDistanceRow( );
						synchronized(this){
							this.row = null;
							synchronized( available ) {
					//			System.out.println("worker "+index+" available");
								available.add( this );
							}
							synchronized( metric ) {
					//			System.out.println("worker "+index+" notified metric");
								metric.notify();
							}
							if(stop){
								return;
							}
							try {
					//			System.out.println("worker "+index+" waits");
								wait();
							} catch ( InterruptedException e ) {
								e.printStackTrace();
							}
						}
					} catch ( Exception e ) {
						e.printStackTrace();
						this.stop = true;
						throw new RuntimeException();
					}
				}
			}
		}
		
		private void fillDistanceRow() throws Exception{
			for(int j=0;j<row.length;j++){
				synchronized(objects[j]){
					row[j] = metric.getDistance( prof1, prof1Rc, objects[j], length );
				}
			}
		}
		
	}
}
