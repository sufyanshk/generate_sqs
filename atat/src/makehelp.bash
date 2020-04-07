#!/bin/bash

echo const char \*helpstring=\"\"
sed 's/^/\"/g' | sed 's/$/\\n\"/g'
echo \;
