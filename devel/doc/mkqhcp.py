#!/usr/bin/python3

## mkqhcp.py
## Version 1.0.4

## Copyright 2022-2023 John Donoghue
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

import sys
import os
import re

def process(name):
  with open(name + ".qhcp", 'wt') as f:
    f.write ('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write ('<QHelpCollectionProject version="1.0">\n')
    f.write ('  <docFiles>\n')
    f.write ('    <generate>\n')
    f.write ('      <file>\n')
    f.write ('        <input>{0}.qhp</input>\n'.format(name))
    f.write ('        <output>{0}.qch</output>\n'.format(name))
    f.write ('      </file>\n')
    f.write ('    </generate>\n')
    f.write ('    <register>\n')
    f.write ('      <file>{0}.qch</file>\n'.format(name))
    f.write ('    </register>\n')
    f.write ('  </docFiles>\n')
    f.write ('</QHelpCollectionProject>\n')

  title = name
  pat_match = re.compile(r".*<title>(?P<title>[^<]+)</title>.*")
  with open(name + ".html", 'rt') as fin:
    # find html
    for line in fin:
      line = line.strip()
      e = pat_match.match(line)
      if e:
          title = e.group("title")
          break

  # section
  h2_match = re.compile(r'.*<h2 class="chapter"[^>]*>(?P<title>[^<]+)</h2>.*')
  # appendix
  h2a_match = re.compile(r'.*<h2 class="appendix"[^>]*>(?P<title>[^<]+)</h2>.*')
  # index
  h2i_match = re.compile(r'.*<h2 class="unnumbered"[^>]*>(?P<title>[^<]+)</h2>.*')

  h3_match = re.compile(r'.*<h3 class="section"[^>]*>(?P<title>[^<]+)</h3>.*')
  h4_match = re.compile(r'.*<h4 class="subsection"[^>]*>(?P<title>[^<]+)</h4>.*')
  tag_match1 = re.compile(r'.*<span id="(?P<tag>[^"]+)"[^>]*></span>.*')
  #tag_match2 = re.compile(r'.*<div class="[sub]*section" id="(?P<tag>[^"]+)"[^>]*>.*')
  tag_match2 = re.compile(r'.*<div class="[sub]*section[^"]*" id="(?P<tag>[^"]+)"[^>]*>.*')
  tag_match3 = re.compile(r'.*<div class="chapter-level-extent" id="(?P<tag>[^"]+)"[^>]*>.*')
  tag_match4 = re.compile(r'.*<div class="appendix-level-extent" id="(?P<tag>[^"]+)"[^>]*>.*')
  tag_match5 = re.compile(r'.*<div class="unnumbered-level-extent" id="(?P<tag>[^"]+)"[^>]*>.*')
  index_match = re.compile(r'.*<h4 class="subsection"[^>]*>[\d\.\s]*(?P<name>[^<]+)</h4>.*')

  tag = "top"
  has_h2 = False 
  has_h3 = False 

  #pat_match = re.compile(r'.*<span id="(?P<tag>[^"])"></span>(?P<title>[.]+)$')
  with open(name + ".html", 'rt') as fin:
    with open(name + ".qhp", 'wt') as f:
      f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
      f.write('<QtHelpProject version="1.0">\n')
      f.write('  <namespace>octave.community.{}</namespace>\n'.format(name))
      f.write('  <virtualFolder>doc</virtualFolder>\n')
      f.write('  <filterSection>\n')
      f.write('    <toc>\n')
      f.write('      <section title="{} Manual" ref="{}.html">\n'.format(title, name))
      # chapters here
      for line in fin:
          line = line.strip()
          e = tag_match1.match(line)
          if not e:
              e = tag_match2.match(line)
          if not e:
              e = tag_match3.match(line)
          if not e:
              e = tag_match4.match(line)
          if not e:
              e = tag_match5.match(line)
          if e:
              tag = e.group("tag")

          e = h2_match.match(line)
          if not e:
              e = h2a_match.match(line)
          if not e:
              e = h2i_match.match(line)
          if e:
              if has_h3:
                  f.write('          </section>\n')
                  has_h3 = False
              if has_h2:
                  f.write('        </section>\n')
              has_h2 = True
              f.write('        <section title="{}" ref="{}.html#{}">\n'.format(e.group("title"), name, tag))

          e = h3_match.match(line)
          if e:
              if has_h3:
                  f.write('          </section>\n')
              has_h3 = True

              f.write('          <section title="{}" ref="{}.html#{}">\n'.format(e.group("title"), name, tag))

          e = h4_match.match(line)
          if e:
              f.write('            <section title="{}" ref="{}.html#{}"></section>\n'.format(e.group("title"), name, tag))


      if has_h3:
        f.write('          </section>\n')
      if has_h2:
        f.write('        </section>\n')

      f.write('      </section>\n')
      f.write('    </toc>\n')
      f.write('    <keywords>\n')
      
      fin.seek(0)
      for line in fin:
          line = line.strip()
          e = tag_match1.match(line)
          if not e:
              e = tag_match2.match(line)
          if e:
              tag = e.group("tag")

          e = index_match.match(line)
          if e:
              f.write('      <keyword name="{}" ref="{}.html#{}"></keyword>\n'.format(e.group("name"), name, tag))

      f.write('    </keywords>\n')
      f.write('    <files>\n')
      f.write('      <file>{}.html</file>\n'.format(name))
      f.write('      <file>{}.css</file>\n'.format(name))
      f.write('    </files>\n')
      f.write('  </filterSection>\n')
      f.write('</QtHelpProject>\n')


def show_usage():
  print (sys.argv[0], "projname")

if __name__ == "__main__":
  if len(sys.argv) > 1:
    status = process(sys.argv[1])
    sys.exit(status)
  else:
    show_usage()
