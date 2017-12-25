#!/usr/bin/env node

matrix = require("./matrix.js");
var argv = require('yargs').argv;

var n = argv.n;
console.log(n);

a = matrix.rand_matrix(n,n);
b = matrix.rand_matrix(n,n);
c = matrix.create_matrix(n,n);

/* Just want to time this part */
console.time('javascript_time');
matrix.mult(c,a,b);
console.timeEnd('javascript_time');

