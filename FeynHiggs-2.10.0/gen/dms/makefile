TAGS = SM Hp Cha Neu Glu

all: f/deltaMs.h

f/deltaMs.h: f/%.h: $(foreach TAG,$(TAGS),m/%-$(TAG).F)
	mkdir -p f
	cat $^ > $@

m/deltaMs-%.F: dmq
	./deltaMs $(@:m/deltaMs-%.F=%)

clean:
	rm -fr m f

.PRECIOUS: %.F

