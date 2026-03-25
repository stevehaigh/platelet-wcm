# Container image #3: wcm-webapp.
# Extends wcm-code with the Dash web UI as the entrypoint.
#
# To build locally from the wcEcoli/ project root directory:
#
#     > docker build -f cloud/docker/webapp/Dockerfile \
#         --build-arg from=${USER}-wcm-code \
#         -t ${USER}-wcm-webapp .
#
# To run locally:
#
#     > docker run --rm -p 8050:8050 -v $(pwd)/out:/wcEcoli/out ${USER}-wcm-webapp
#
# Then open http://localhost:8050/

ARG from=steve-wcm-code:latest
FROM ${from}

# Copy webapp source files (not present in the base wcm-code image)
COPY wholecell/webapp/ /wcEcoli/wholecell/webapp/
COPY runscripts/manual/webapp.py /wcEcoli/runscripts/manual/webapp.py

# The webapp uses local subprocess for simulations (no docker-in-docker).
ENV WCECOLI_WEBAPP_MODE=container

EXPOSE 8050

CMD ["python", "runscripts/manual/webapp.py", "--port", "8050"]
