This is the LaTeX template for KMD Master's thesis (2015 edition).

1. ABOUT LATEX

LaTeX is a typesetting and document publishing system (NOT a word processor). 
Basically, you write the document contents in plain text, using a kind of 
markup language, then generate the document (PDF) using LaTeX and the provided
template.

The idea is that you can focus on writing the CONTENTS, and not care (much) 
about the formatting, typography, page layout, citation style, etc.  Another
advantage is that you can divide your document into many different files (for
example, a separate file for each chapter and appendix), and edit only the
part which you want to work on without loading the entire document at once.

A functional LaTeX system is actually a large collection of programs, macros, 
templates, and support files which work together.  So you need to install a 
LaTeX _distribution_ which bundles it all together and sets it up on your 
system.  We recommend using W32TeX for Windows, MacTeX for MacOS, or TeXLive 
for Linux.

A useful introductory reference for LaTeX: http://en.wikibooks.org/wiki/LaTeX


2. INSTALLING LATEX

For Windows (XP/Vista/7/8):

  Get W32TeX using "TeX Installer 3":
  http://www.math.sci.hokudai.ac.jp/%7Eabenori/soft/bin/abtexinst_0_84r6.zip

  Unfortunately the installer is in Japanese only.  Just accept all the default
  settings, EXCEPT that on the second wizard page, the second checkbox 
  ("ネットワークを使わない．（ローカルのファイルからインストールします．）") should be UN-checked
  (unless you are installing from a previously-downloaded distribution).

For MacOS X:

  Get MacTeX 2014 from: http://www2.kmd.keio.ac.jp/~kato/mactex/
  Follow the instructions on the page.

For Linux/BSD:

  Use your distribution's package manager to install TeXLive.


3. USING THE TEMPLATE FILES

Unzip the files into a new directory.  The files ending in .sty and .bst are
template files containing LaTeX code; do not modify them, but simply place 
them in the same directory as your document files.  The files ending in .tex 
and .bib are the sample document files - they are plain text files which can 
be edited with any text editor.  Use them as a basis for making your own 
document:
  document.tex     OR  document-en.tex       The main document (start here)
  intro.tex        OR  intro-en.tex          Contents of sample chapter 1 
  relatedworks.tex OR  relatedworks-en.tex   Contents of sample chapter 2
  appendix1.tex    OR  appendix1-en.tex      Contents of sample appendix A
  newchapter.tex   OR  newchapter-en.tex     Empty template for a new chapter
  document.bib     OR  document-en.bib       List of references (bibliography)
  figures\*.jpg                              Images used in the document

You can choose to start appendices with \section commands (as in the sample 
appendix here), or with \chapter commands. However, you must use the same 
method for all your appendices. (Using \chapter for each appendix will cause 
it to start on a new page, which might be suitable if you have very long 
appendices, or a large number of them.)

The file 'useful_examples.txt' contains various useful examples of techniques 
for lists, tables, figures, and so on that you can copy into your own thesis 
and modify.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Note For Windows/W32TeX Users:

The W32TeX distribution on Windows does not include the 'cite', 'subfigure' and
'tabu' packages by default.  If you don't need these packages, you can comment
them out of the main .tex file.  Otherwise, extract the cite.sty, subfigure.sty
and tabu.sty from w32tex_extra.zip and place them in the same directory as the 
other template and document files.  

w32tex_extra also contains the file 'cid-x.map', which is a modified version
of the dvipdfmx font map file for Asian fonts.  You are recommended to copy
this file to your \w32tex\share\texmf-dist\fonts\map\dvipdfmx\base directory
(replacing the existing file). This changes the default Japanese fonts to 
default PostScript fonts, which allows bold headings in Japanese (and also 
reduces the PDF file size somewhat).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


USING OVERLEAF

Overleaf (a.k.a. WriteLaTeX) is an on-line web-based editor for LaTeX:
https://www.overleaf.com/

To use it, you have to create an account on their web site (or log in via 
Google), but a basic account is free, and includes all the features you should 
need for writing your thesis.  (The basic account allows a maximum of 60 files 
per project, but unless you have a very large number of images, that should be 
plenty.}

To start a new project based on this thesis template, simply choose the upload
(arrow) on the dashboard after signing in, and choose 'Upload Zip'.  Then drag
and drop the "KMD_New_Thesis_Files_(...).zip" file to upload it.  You will see
a warning that w32tex_extra.zip could not be uploaded -- this is fine, so just
ignore that message. A new KMD thesis project will then be created.

If you're using the Japanese template, you must go into the settings (gear icon
on the top right), and under Project Settings, Advanced Build Options, change 
the 'LaTeX Engine' selection to "LaTeX dvipdf".  Then choose "Save Project 
Settings".  Until you do this, you will get an error when compiling.
