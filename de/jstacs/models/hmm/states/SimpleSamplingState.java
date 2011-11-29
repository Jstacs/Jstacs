/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.hmm.states;

import java.io.IOException;

import de.jstacs.models.hmm.states.emissions.SamplingEmission;

/**
 * This class implements a state that can be used for a HMM that obtains its parameters from sampling.
 * 
 * @author Jens Keilwagen
 */
public class SimpleSamplingState extends SimpleState implements SamplingState {

	/**
	 * This constructor creates a state that can be used in a HMM that obtains its parameters from sampling.
	 * 
	 * @param e the emission
	 * @param name the name of the emission (e.g. for the Graphviz representation)
	 * @param forward a switch whether to use the forward or the reverse strand
	 */
	public SimpleSamplingState( SamplingEmission e, String name, boolean forward ) {
		super( e, name, forward );
	}
	
	public void drawParametersFromStatistic() throws Exception {
		((SamplingEmission)e).drawParametersFromStatistic();		
	}

	public void extendSampling( int sampling, boolean append ) throws IOException {
		((SamplingEmission)e).extendSampling( sampling, append );
	}

	public void initForSampling( int starts ) throws IOException {
		((SamplingEmission)e).initForSampling( starts );		
	}

	public boolean isInSamplingMode() {
		return ((SamplingEmission)e).isInSamplingMode();
	}

	public boolean parseNextParameterSet() {
		return ((SamplingEmission)e).parseNextParameterSet();
	}

	public boolean parseParameterSet( int sampling, int n ) throws Exception {
		return ((SamplingEmission)e).parseParameterSet( sampling, n );
	}

	public void samplingStopped() throws IOException {
		((SamplingEmission)e).samplingStopped();		
	}

    public void acceptParameters() throws IOException {
        ((SamplingEmission)e).acceptParameters();
    }
    
    public double getLogGammaScoreForCurrentStatistic() {
		return ((SamplingEmission)e).getLogGammaScoreFromStatistic();
	}
    
	public double getLogPosteriorFromStatistic() {
		return ((SamplingEmission)e).getLogPosteriorFromStatistic();
	}
}
