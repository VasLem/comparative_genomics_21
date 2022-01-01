# as adapted from https://github.com/gamcil/synthaser/blob/8299373b997f175f536628720e4c66bbff42b4ab/synthaser/results.py
import ncbi
import logging


LOG = logging.getLogger(__name__)


def search(
    queries: list,
    results_file=None,
    cdsid=None,
    delay=20,
    max_retries=-1,
    database=None,
    **kwargs,
):
    query = '\n'.join([f'>{cnt}\n{q}' for cnt, q in enumerate(queries)])

    if len(query) > 4000:
        raise ValueError("Too many sequences (NCBI limit = 4000)")

    handle = _remote(
            query,
            output=results_file,
            cdsid=cdsid,
            delay=delay,
            max_retries=max_retries,
            database=database,
            **kwargs,
        )

    return handle


def _remote(query, cdsid=None, delay=20, max_retries=-1, output=None, **kwargs):
    """Launch new CD-Search job, poll results and return a faux 'handle' for parsing."""

    ncbi.set_search_params(**kwargs)

    if not cdsid:
        LOG.info("Launching new CD-Search run")
        cdsid = ncbi.launch(query)

    LOG.info("Run ID: %s", cdsid)
    LOG.info("Polling NCBI for results...")
    response = ncbi.retrieve(cdsid, delay=delay, max_retries=max_retries)
    if output:
        LOG.info("Writing CD-Search results table to %s", output)
        with open(output, "w") as out:
            out.write(response.text)

    return response.text.split("\n")
