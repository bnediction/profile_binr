
EXPOSE_PORT = 8888

run: build
	docker run -it --rm --volume $(PWD):/notebook -w /notebook -p $(EXPOSE_PORT):8888 bnediction/profile_binr:main

run-isolated: build
	docker run -it --rm -w /notebook -p $(EXPOSE_PORT):8888 bnediction/profile_binr:main


build:
	docker build -t bnediction/profile_binr:local .


