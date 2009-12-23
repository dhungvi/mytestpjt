/*
 * FactorSequence.java
 *
 * Copyright 2003 Sergio Anibal de Carvalho Junior
 *
 * This file is part of NeoBio.
 *
 * NeoBio is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * NeoBio is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with NeoBio;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Proper attribution of the author as the source of the software would be appreciated.
 *
 * Sergio Anibal de Carvalho Junior		mailto:sergioanibaljr@users.sourceforge.net
 * Department of Computer Science		http://www.dcs.kcl.ac.uk
 * King's College London, UK			http://www.kcl.ac.uk
 *
 * Please visit http://neobio.sourceforge.net
 *
 * This project was supervised by Professor Maxime Crochemore.
 *
 */

package neobio.alignment;

import java.io.Reader;
import java.io.BufferedReader;
import java.io.IOException;

/**
 * This class builds a list of factors of a character sequence as induced by its
 * Lempel-Ziv (LZ78) factorisation. Each factor is enconded as the longest factor
 * previously seen plus one character.
 *
 * <P>The input can come from any source, provided it is encapsulated in a proper
 * <CODE>Reader</CODE> instance. The stream is expected to be ready (i.e. the next
 * <CODE>read</CODE> operation must return the first character of the sequence) and it is
 * not closed when its end is reached, so the client is allowed to reset it and maybe use
 * it for another purpose.</P>
 *
 * <P>Sequences can contain letters only although lines started with the
 * <CODE>COMMENT_CHAR</CODE> character ('>') are regarded as comments and are completely
 * skipped. White spaces (including tabs, line feeds and carriage returns) are also
 * ignored throughout.</P>
 *
 * <P>This class uses a {@linkplain Trie} to keep track of a list of factors. Each node of
 * the trie contains a {@linkplain Factor} of the text. As the sequence is read from the
 * input, the trie is traversed as far as possible. When a leaf node is reached (which
 * means that the longest prefix of the input has been found), two tasks are
 * accomplished:</P>
 *
 * <UL>
 * <LI>a new <CODE>Factor</CODE> is created with the character at the current position of
 * the input and the leaf node's factor;
 * <LI>a new node is added to the trie with the character at the current position of the
 * input;
 * </UL>
 *
 * <P>Each factor also receives a serial number according to the order they are found and
 * a pointer to the next factor (in that order) for fast access. This pointer, together
 * with the factor's ancestor pointer forms a doubly-linked list of factors. The original
 * text can then be reconstructed simply by following the linked list and writing out its
 * factors.</P>
 *
 * <P>As an example, the sequence <CODE>ACTAAACCGCATTAATAATAAAA</CODE> is parsed into the
 * following 12 factors:</P>
 *
 * <CODE><BLOCKQUOTE><PRE>
 * 0  ( , ) = empty
 * 1  (0,A) = A
 * 2  (0,C) = C
 * 3  (0,T) = T
 * 4  (1,A) = AA
 * 5  (1,C) = AC
 * 6  (2,G) = CG
 * 7  (2,A) = CA
 * 8  (3,T) = TT
 * 9  (4,T) = AAT
 * 10 (9,A) = AATA
 * 11 (4,A) = AAA
 *
 * serial # (prefix, new char) = factor text
 * </PRE></BLOCKQUOTE></CODE>
 *
 * <P>This class is used by {@linkplain CrochemoreLandauZivUkelson} algorithm to speed up
 * the classic dynamic programming approach to sequence alignment.</P>
 *
 * @author Sergio A. de Carvalho Jr.
 * @see Factor
 * @see Trie
 * @see CrochemoreLandauZivUkelson
 */
public class FactorSequence
{
	/**
	 * The character used to start a comment line in a sequence file. When this character
	 * is found, the rest of the line is ignored.
	 */
	protected static final char COMMENT_CHAR = '>';

	/**
	 * A pointer to the root factor, the one that starts the list of factors.
	 */
	protected Factor root_factor;

	/**
	 * The numbers of character represented by this sequence.
	 */
	protected int num_chars;

	/**
	 * The numbers of factors generated by the LZ78 parsing of the sequence.
	 */
	protected int num_factors;

	/**
	 * Creates a new instance of a <CODE>FactorSequence</CODE>, loading the sequence data
	 * from the <CODE>Reader</CODE> input stream. A doubly-linked list of factors is built
	 * according to its LZ78 factorisation.
	 *
	 * @param reader source of characters for this sequence
	 * @throws IOException if an I/O exception occurs when reading the input
	 * @throws InvalidSequenceException if the input does not contain a valid sequence
	 */
	public FactorSequence (Reader reader)
		throws IOException, InvalidSequenceException
	{
		BufferedReader	input = new BufferedReader(reader);
		Trie			root_node, current_node, new_node = null;
		Factor			current_factor, last_factor, new_factor;
		int				ch;
		char			c;

		// create root factor and the root node of the trie
		root_factor = new Factor ();
		root_node = new Trie (root_factor);
		num_factors = 1;
		num_chars = 0;

		current_node = root_node;
		last_factor = root_factor;

		// read characters from the input
		while ((ch = input.read()) != -1)
		{
			c = (char) ch;

			if (c == COMMENT_CHAR)
				// it's a comment line: skip it!
				input.readLine();

			// accept letters only
			else if (Character.isLetter(c))
			{
				num_chars++;

				// walk down the trie as far as possible
				new_node = current_node.spellDown(c);

				if (new_node != null)
				{
					current_node = new_node;
				}
				else
				{
					// the longest factor of the input has been found,
					// now create a new factor from the current node's factor
					current_factor = (Factor) current_node.getData();
					new_factor = new Factor (current_factor, num_factors, c);

					// add the new character to the trie as well
					current_node.add (new_factor, c);

					// set up a pointer from the last factor to the new one
					last_factor.setNext (new_factor);
					last_factor = new_factor;

					// restart at the root of the trie
					current_node = root_node;

					num_factors++;
				}
			}

			// anything else, except whitespaces, will throw an exception
			else if (!Character.isWhitespace(c))
				throw new InvalidSequenceException
					("Sequences can contain letters only.");
		}

		// if new_node is not null, the last factor is actually
		// not a new factor but a factor already created
		if (new_node != null)
		{
			// no new node is created, just point the last_factor to an
			// existing one that represents the last characters of the text
			last_factor.setNext((Factor) new_node.getData());

			num_factors++;
		}

		// check if read anything useful!
		if (num_factors <= 1)
			throw new InvalidSequenceException ("Empty sequence.");
	}

	/**
	 * Returns the root factor, the one that starts the list of factors.
	 *
	 * @return root factor
	 */
	public Factor getRootFactor ()
	{
		return root_factor;
	}

	/**
	 * Returns the number of factors produced by the LZ78 parsing of the text.
	 *
	 * @return number of factors
	 */
	public int numFactors()
	{
		return num_factors;
	}

	/**
	 * Returns the number of characters of the original sequence.
	 *
	 * @return number of characters of the original sequence
	 */
	public int numChars ()
	{
		return num_chars;
	}

	/**
	 * Reconstructs the sequence from the list of factors induced by the LZ78 parsing of
	 * the text.
	 *
	 * @return the original sequence
	 */
	public String toString ()
	{
		StringBuffer	buf = new StringBuffer();
		Factor			node;

		node = root_factor.getNext();

		for (int i = 1; i < numFactors(); i++)
		{
			buf.append(node);

			node = node.getNext();
		}

		return buf.toString();
	}

	/**
	 * Returns a string representation of the actual list of factors produced by the LZ78
	 * parsing of the text. Each factor is printed out in a separate line, in the order
	 * they appear in the text, with its serial number, its ancestor's serial number, its
	 * new character, length and a string representation of the factor itself.
	 *
	 * @return a string representation of the list of factors
	 */
	public String printFactors ()
	{
		StringBuffer	buf = new StringBuffer();
		Factor			factor;

		factor = root_factor.getNext();

		for (int i = 1; i < numFactors(); i++)
		{
			buf.append (factor.getSerialNumber() + "\t<");
			buf.append (factor.getAncestor().getSerialNumber() + " ,\t");
			buf.append (factor.getNewChar() + ">\t");
			buf.append (factor.length() + "\t" + factor + "\n");

			factor = factor.getNext();
		}

		buf.append(numFactors() + " factors\n");

		return buf.toString();
	}
}
