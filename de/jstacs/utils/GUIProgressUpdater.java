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

package de.jstacs.utils;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JProgressBar;

/**
 * This class implements a {@link ProgressUpdater} with a GUI.
 * 
 * @author Jens Keilwagen
 */
public class GUIProgressUpdater implements ProgressUpdater, ActionListener {

	private JProgressBar progress;

	private JFrame frame;

	private JButton b;

	private boolean canceled;

	/**
	 * This is the constructor for a {@link GUIProgressUpdater}. It enables the
	 * programmer to add a cancel button.
	 * 
	 * @param cancelButton
	 *            <code>true</code> if the cancel button shall be added
	 */
	public GUIProgressUpdater( boolean cancelButton ) {
		progress = new JProgressBar();
		progress.setIndeterminate( false );
		progress.setStringPainted( true );
		progress.setMinimum( 0 );

		if( cancelButton ) {
			b = new JButton( "cancel" );
			b.addActionListener( this );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setMax(int)
	 */
	public void setMax( int max ) {
		progress.setMaximum( max );
		frame = new JFrame();
		frame.setDefaultCloseOperation( JFrame.DO_NOTHING_ON_CLOSE );
		frame.setResizable( false );
		frame.setLayout( new FlowLayout() );
		frame.add( progress );
		if( b != null ) {
			frame.add( b );
		}
		frame.pack();
		frame.setVisible( true );
		canceled = false;
		setValue( 0 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#setValue(int)
	 */
	public void setValue( int value ) {
		progress.setValue( value );
		frame.setTitle( progress.getString() + " done" );
		if( progress.getMaximum() == value ) {
			try {
				Thread.sleep( 500 );
			} catch ( InterruptedException e ) {}
			frame.dispose();
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		frame.dispose();
		super.finalize();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.utils.ProgressUpdater#isCancelled()
	 */
	public boolean isCancelled() {
		return canceled;
	}

	/* (non-Javadoc)
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed( ActionEvent e ) {
		canceled = true;
		frame.setVisible( false );
		frame.dispose();
	}
}
