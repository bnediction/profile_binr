function updater
  Rscript tidy_source.R
  poetry run python updateR.py
end
