#!/bin/bash
set -Ceu

echo 'Do you really remove the followig?'
ls outputs/
read -p '' 

rm outputs/*
