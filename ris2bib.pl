#!/usr/bin/perl
#
# Convert RIS data into BibTeX
#
# (C) Copyright 2005-2006 Diomidis Spinellis
#
# Permission to use, copy, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies and that
# both that copyright notice and this permission notice appear in
# supporting documentation.
#
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#


# Read RIS records from the standard input producing BibTeX records
# on the standard output.

while (<>) {
	chop;
	if (/^TY\s*\-\s*(.*)/) {
		$type = $1;
		flush_record();
		$rtype = "InProceedings" if ($type eq 'CONF');
		$rtype = "Article" if ($type eq 'JOUR');
		$rtype = "Book" if ($type eq 'BOOK');
		undef $author;
		undef $title;
		undef $start_page;
		undef $end_page;
		undef $year;
		undef $doi;
		undef $url;
		undef $volume;
		undef $number;
		undef $pages;
		undef $journal;
		undef $booktitle;
	}
	if (/^AU\s*\-\s*(.*)/) {
		if (defined($author)) {
			$author .= " and $1";
		} else {
			$author = $1;
		}
	}
	$title = $1 if (/^T1\s*\-\s*(.*)/);
	if (/^JF\s*\-\s*(.*)/ && $rtype eq 'Book') {
		$rtype = 'InCollection';
		$booktitle = $1;
	}
	$start_page = $1 if (/^SP\s*\-\s*(.*)/);
	$end_page = $1 if (/^EP\s*\-\s*(.*)/);
	$year = $1 if (/^PY\s*\-\s*(\d+)/);
	# Scopus DOI
	$doi = $1 if (/^N1\s*\-\s*doi\:\s*(.+)/i);
	$url = $1 if (/^UR\s*\-\s*(.+)/);
	# SprignerLink DOI
	$doi = $1 if (/^M3\s*\-\s*(.+)/i);
	if (/^\s*J[FO]\s*\-\s*(.*)/) {
		if ($rtype eq 'InProceedings') {
			$booktitle = $1;
		} elsif ($rtype eq 'Article') {
			$journal = $1;
		}
	}
	$volume = $1 if (/^VL\s*\-\s*(.*)/);
	$number = $1 if (/^IS\s*\-\s*(.*)/);
}

flush_record();

sub
flush_record
{
	return unless defined($rtype);
	$n++;
	if ($start_page) {
		$pages = $start_page;
		$pages .= "--$end_page" if ($end_page);
	}
	print "\@$rtype\{EXPORT$n,\n";
	print "\tAuthor={$author},\n" if ($author);
	print "\tTitle={$title},\n" if ($title);
	print "\tBooktitle={$booktitle},\n" if ($booktitle);
	print "\tJournal={$journal},\n" if ($journal);
	print "\tPages={$pages},\n" if ($pages);
	print "\tYear={$year},\n" if ($year);
	print "\tDOI={$doi},\n" if ($doi);
	print "\tVolume={$volume},\n" if ($volume);
	print "\tNumber={$number},\n" if ($number);
	print "\tURL={$url},\n" if ($url);
	print "\tXExportBy={ris2bib}\n}\n\n";
}

