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

package de.jstacs.classifier.differentiableSequenceScoreBased;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.data.DataSet;

/**
 * This class enables the user to exploit all CPUs of an computer by using threads.
 * Each thread computes the function for a part of the data.
 * The number of compute threads can be determined in the constructor.
 * 
 * <br>
 * <br>
 * 
 * The main goal is on the one hand to hide the complete infrastructure that is used for the multi-threading from
 * the programmer and on the other hand to provide an simple interface in form of few abstract methods.
 * 
 * <br>
 * <br>
 * 
 * It is very important for this class that the clone() method (of the used items) works correctly, since each thread works on its own clones. 
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractMultiThreadedOptimizableFunction extends AbstractOptimizableFunction
{
	/**
	 * This method returns the number of available processors. So if no other job is running this is the recommended number of threads.
	 * 
	 * @return the number of available processors
	 */
	public static final int getNumberOfAvailableProcessors()
	{
		return Runtime.getRuntime().availableProcessors();
	}

	private Worker[] worker;
	
	/**
	 * This is a pointer for the current parameters.
	 */
	protected double[] params;

	/**
	 * The constructor for an multi-threaded instance.
	 * 
	 * @param threads the number of threads used for evaluating the function and determining the gradient of the function
	 * @param data the array of {@link DataSet}s containing the data that is needed to evaluate the function
	 * @param weights the weights for each {@link de.jstacs.data.sequences.Sequence} in each {@link DataSet} of  <code>data</code>
	 * @param norm
	 *            the switch for using the normalization (division by the number
	 *            of sequences)
	 * @param freeParams
	 *            the switch for using only the free parameters
	 * 
	 * @throws IllegalArgumentException
	 *             if the number of threads is not positive, the {@link DataSet}s contain too less {@link de.jstacs.data.sequences.Sequence} for the number of threads, the number of classes or the dimension of the weights is not correct
	 */
	public AbstractMultiThreadedOptimizableFunction( int threads, DataSet[] data, double[][] weights, boolean norm, boolean freeParams ) throws IllegalArgumentException
	{
		super( data, weights, norm, freeParams );
		if( threads < 1 )
		{
			throw new IllegalArgumentException( "The number of threads has to be positive." );
		}
		worker = new Worker[threads];
		prepareThreads();
	}
	
	
	
	public void setDataAndWeights( DataSet[] data, double[][] weights ) throws IllegalArgumentException {
		super.setDataAndWeights(data, weights);
		if( worker != null ) {
			prepareThreads();
		}
	}
	
	//assigns parts of the data to the threads
	private void prepareThreads() {
		int i = 0, anz = 0;
		for( ; i < data.length; i++ )
		{
			anz += data[i].getNumberOfElements();
		}
		anz = (int) Math.ceil( anz / (double) worker.length );
		int startClass, endClass = 0, startSeq, endSeq = 0, c; 
		for( i = 0; i < worker.length; i++ )
		{
			startSeq = endSeq;
			startClass = endClass;
			if( i == worker.length-1 )
			{
				endClass = data.length-1;
				endSeq = data[endClass].getNumberOfElements();
			}
			else
			{
				c = anz;
				while( endClass < data.length && data[endClass].getNumberOfElements()- endSeq < c )
				{
					c -= (data[endClass].getNumberOfElements()- endSeq);
					endSeq = 0;
					endClass++;
				}
				endSeq += c;
				if( endClass >= data.length || startClass == endClass && startSeq == endSeq ) {
					throw new IllegalArgumentException( "There are less sequence than threads used for the optimization. This seems to be unlikely. Please check your data or reduce the number of threads." );
				}
			}
			//System.out.println("split " + j + ": " + startClass + " " + startSeq + "\t" + endClass + " " + endSeq );
			if( worker[i] != null ) {
				if( worker[i].isWaiting() ) {
					worker[i].setIndices(startClass, startSeq, endClass, endSeq);
				} else {
					stopThreads();
					throw new RuntimeException();
				}
			} else {
				worker[i] = new Worker( i, startClass, startSeq, endClass, endSeq );
				worker[i].start();
			}
		}
	}

	public final double[] evaluateGradientOfFunction( double[] x ) throws DimensionException, EvaluationException
	{
		setParams( x );
		waitUntilWorkersFinished( WorkerTask.EVALUATE_GRADIENT );
		return joinGradients();
	}
	
	/**
	 * This method evaluates the gradient of the function for a part of the data.
	 * 
	 * @param index the index of the part
	 * @param startClass the index of the start class
	 * @param startSeq the index of the start sequence
	 * @param endClass the index of the end class (inclusive)
	 * @param endSeq the index of the end sequence (exclusive)
	 */
	protected abstract void evaluateGradientOfFunction( int index, int startClass, int startSeq, int endClass, int endSeq );
	
	/**
	 * This method joins the gradients of each part that have been computed using {@link AbstractMultiThreadedOptimizableFunction#evaluateGradientOfFunction(int, int, int, int, int)}.
	 * 
	 * @return the gradient
	 * 
	 * @throws EvaluationException if the gradient could not be evaluated properly 
	 */
	protected abstract double[] joinGradients() throws EvaluationException;
	
	public final double evaluateFunction( double[] x ) throws DimensionException, EvaluationException
	{
		setParams( x );
		waitUntilWorkersFinished( WorkerTask.EVALUATE );
		return joinFunction();
	}
	
	/**
	 * This method evaluates the function for a part of the data.
	 * 
	 * @param index the index of the part
	 * @param startClass the index of the start class
	 * @param startSeq the index of the start sequence
	 * @param endClass the index of the end class (inclusive)
	 * @param endSeq the index of the end sequence (exclusive)
	 * 
	 * @throws EvaluationException if the gradient could not be evaluated properly
	 */
	protected abstract void evaluateFunction( int index, int startClass, int startSeq, int endClass, int endSeq ) throws EvaluationException;
	
	/**
	 * This method joins the partial results that have been computed using {@link AbstractMultiThreadedOptimizableFunction#evaluateFunction(int, int, int, int, int)}.

	 * @return the value of the function
	 * 
	 * @throws EvaluationException if the gradient could not be evaluated properly
	 * @throws DimensionException if the parameters could not be set
	 */
	protected abstract double joinFunction() throws EvaluationException, DimensionException;
	
	public final void setParams( double[] params ) throws DimensionException
	{
		if( this.params == null || this.params.length != params.length ) {
			this.params = params.clone();
		} else {
			System.arraycopy( params, 0, this.params, 0, params.length );
		}
		setThreadIndependentParameters();
		waitUntilWorkersFinished( WorkerTask.SET_PARAMETERS );
	}
	
	/**
	 * This method allows to set thread independent parameters.
	 * It also allows to check some constraints as for instance the dimension of the parameters.
	 * 
	 * @throws DimensionException if the dimension of the parameters is wrong
	 * 
	 * @see #setParams(double[])
	 */
	protected abstract void setThreadIndependentParameters() throws DimensionException;
	
	/**
	 * This method sets the parameters for thread <code>index</code>
	 * 
	 * @param index the index of the thread
	 * 
	 * @throws DimensionException if the parameters could not be set
	 * 
	 * @see AbstractMultiThreadedOptimizableFunction#params
	 */
	protected abstract void setParams( int index ) throws DimensionException;

	/**
	 * This method waits until all worker (threads) have done the specified {@link WorkerTask}.
	 * 
	 * @param wt the task
	 */
	private synchronized void waitUntilWorkersFinished( WorkerTask wt )
	{
		for( int t = 0; t < worker.length; t++ )
		{
			worker[t].setTask( wt );
		}
		
		boolean exception = false;
		int t = -1, i;
		while( true )
		{
			i = 0;
			while( i < worker.length && worker[i].isWaiting() ){
				if( worker[i].exception ) {
					exception = true;
					t = i;
				}
				i++;
			}
			if( i == worker.length ){
				if( exception ) {
					for( i = 0; i < worker.length; i++ ) {
						worker[i].interrupt();
					}
					stopThreads();
					throw new RuntimeException( "Terminate program, since at least thread " + t + " throws an exception." );
				} else {
					break;
				}
			}else{
				try{
					wait();
				}
				catch( InterruptedException e ) { }
			}
		}
	}
	
	/**
	 * This method can and should be used to stop all threads if they are not needed any longer.
	 */
	public final void stopThreads()
	{
		if( worker.length > 1 )
		{
			for( int i = 0; i < worker.length; i++ )
			{
				worker[i].setTask( WorkerTask.STOP );
			}
			worker = null;
		}
	}

	/**
	 * Returns the number of used threads for evaluating the function and for determining the gradient of the function.
	 * 
	 * @return the number of used threads for evaluating the function and for determining the gradient of the function
	 */
	public final int getNumberOfThreads() {
		return worker.length;
	}
	
	/**
	 * This enum defines the task that a worker (thread) has to do.
	 * 
	 * @author Jens Keilwagen
	 *
	 * @see Worker
	 */
	private static enum WorkerTask {
		/**
		 * Indicates that the worker should stop.
		 */
		STOP,
		/**
		 * Indicates that the worker should wait.
		 */
		WAIT,
		/**
		 * Indicates that the worker should set new parameters.
		 */
		SET_PARAMETERS,
		/**
		 * Indicates that the worker should evaluate the function.
		 */
		EVALUATE,
		/**
		 * Indicates that the worker should evaluate the gradient of the function.
		 */
		EVALUATE_GRADIENT;
	}
	
	/**
	 * This class is used to split the computation. Each worker (thread) has only to compute a part of the whole computation.
	 * 
	 * @author Jens Keilwagen
	 */
	private class Worker extends Thread
	{		
		private WorkerTask task;
		private int index, startClass, startSeq, endClass, endSeq;
		private boolean exception;
		
		/**
		 * This constructor creates a {@link Worker} that is used to evaluate the function and its gradient.
		 * The instance only computes the values for a part of the data that is given by specific indices.
		 * 
		 * @param index the index of the worker
		 * @param startClass the start index of the classes 
		 * @param startSeq the start index of the sequences
		 * @param endClass the end index of the classes
		 * @param endSeq the end index of the sequences
		 */
		public Worker( int index, int startClass, int startSeq, int endClass, int endSeq )
		{
			super( "worker thread " + index );
			setDaemon( true );
			this.index = index;
			setIndices(startClass, startSeq, endClass, endSeq);
		}
		
		private void setIndices( int startClass, int startSeq, int endClass, int endSeq ) {
			this.startClass = startClass;
			this.startSeq = startSeq;
			this.endClass = endClass;
			this.endSeq = endSeq;
			task = WorkerTask.WAIT;
		}

		public synchronized void run()
		{
			exception = false;
			while( task != WorkerTask.STOP ){
				if( task != WorkerTask.WAIT ){
					try{
						switch( task )
						{
							case SET_PARAMETERS:
								setParams( index );
								break;
							case EVALUATE:
								evaluateFunction( index, startClass, startSeq, endClass, endSeq );
								break;
							case EVALUATE_GRADIENT:
								evaluateGradientOfFunction( index, startClass, startSeq, endClass, endSeq );
								break;
						}
					}catch( Exception e ){
						exception = true;
						e.printStackTrace();
					}
					synchronized ( AbstractMultiThreadedOptimizableFunction.this ) {
						task = WorkerTask.WAIT;
						AbstractMultiThreadedOptimizableFunction.this.notify();
					}
				} else {
					try {
						wait();
					} catch( InterruptedException e ){}
				}
			}
		}
		
		/**
		 * This method can be used to tell the worker which task he has to perform.
		 * 
		 * @param task the task
		 */
		public synchronized void setTask( WorkerTask task ){
			this.task = task;
			notify();
		}
		
		/**
		 * This method answers the question whether the worker waits for new tasks.
		 *  
		 * @return <code>true</code> if the thread is waiting.
		 * 
		 * @see WorkerTask#WAIT
		 */
		public boolean isWaiting(){
			return task == WorkerTask.WAIT;
		}
	}
}
