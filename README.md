# Preliminary AI investigations into massive modes in M-theory: a method write-up
## Chris Blair

## 1. Summary

This is a method write-up of my prototype investigations with AI into the massive modes of type IIA string theory and their hidden M-theory structure. This has involved trialling different frontier AI models, independently or in conjunction, in order to:
* generate code and data
* formulate analytical arguments using mathematics and physics
* write status reports in markdown and LaTeX
* perform literature searches
* conduct audits of the above
* suggest improvements to the workflow e.g. specific agentic instructions, helper scripts

The starting point was the following prompt given to Gemini 3.7 running in the Antigravity IDE/agent platform:

> The goal here is to study the paper https://arxiv.org/abs/hep-th/9503228v1 and test whether the proposals in III.B and III.C work to all levels

The paper cited is due to Bars, and conjectures the form of extra M-theory structure which would complete the massive string states (of the type IIA string) into consistent M-theoretic states. Bars checked this structure works up to level 5, out of infinity. 

The principal success of my preliminary investigations has been to push the check up to level 40 (with Claude Code) and in doing so identify ambiguities which appear beyond level 5. These higher-level checks have been carried out independently by different models, and all overlaps agree. Some analytic avenues to an all-levels proof have been explored, but the status of these remains unresolved and largely unverified.

## 2. Setup and workflow evolution

| Phase | Dates (2026) | Models | Repository | Main output |
|--|--|--|--|--|
| 1 | late Aug | Gemini 3.7 (Antigravity), ChatGPT (cloud) |  [massive-modes](https://github.com/blaircda/massive-modes) | Bars conjecture to level 10 in plain Python; 19-page LaTeX report |
| 2 | early Sep- mid Oct | Gemini 3.7 (Antigravity), GPT-5.6 Terra/Sol (Codex) |  [wiser-massive-modes](https://github.com/blaircda/massive-modes) | Bars conjecture to level 20 in Python+SageMath, initial computational structure results, copius reports on partial analytic results including dead ends |
| 3 | 30 Sep – mid Oct | Claude Opus 5.5 (Claude Code) | [simpler-massive-modes](https://github.com/blaircda/simpler-massive-modes) | Independent Python+SageMath pipeline to level 40, cleaner data setup |

### Phase 1. Awkward but promising 

- **Setup:** Gemini agents producing Python code running locally, with ChatGPT (in browser) used to analyse results and brainstorm methods
- **Initial work:**
    - With a few false starts, Gemini produced Python code from scratch capable of handling the needed representation theoretic methodology
    - ChatGPT proved very capable of offering guidance and feedback on the code/data produced, catching and correcting various errors
    - The AI succeed in producing a check of the Bars conjecture up to level 10, which I verified independently using SageMath

I viewed this initial phase as a successful experiment but decided that additional work should build directly on SageMath, rather than custom mathematical code, and have a more structured format.
It was also abundantly clear that I was using an awkward workflow with the seemingly more capable AI model restricted to chat sessions in browser.
I decided to go for a clean restart.

### Phase 2. More structure, more progress 

- **Setup:** new project repository, with more refined initial guidance:
    - I added a protected baseline of human-verified scripts containing the level 10 results in SageMath, declared as read-only to agents
    - I wrote an agent skill file for carrying out literature searches using arXiv/Inspire
    - I replaced browser ChatGPT sessions with GPT 5.6 Terra/Sol in Codex
    - I enforced a structured workflow with a central STATUS.md report file and subdirectories for independent research directions, each containing their own local STATUS.md report, with instructions that each agent run finish by updating these files
- **Initial work :**
    - Following a literature search (combined human and AI), I began efforts to understand possible analytic structures which could be exploited for an all-levels proof 
    - Tasking Gemini with more advanced work led to erratic output. In status reports it over-confidently reported successes and definite conclusions based only on tenuous evidence. When tasked with rederiving known results from an alternative starting point, its code led to an infinite loop. To circumvent this, it imported the known results, checked they equalled themselves, and then reported success.
- **Audit and switch to GPT model:**
    - An audit by the GPT model highlighted the errors and unfounded claims in Gemini's codes and reports - these were archived and clearly labelled with correction notices
    - I modified the global agent instructions, adding key rules mandating that i) baseline results are kept separate from new derivations, ii) failures must be reported as obstacles rather than viewed as something to circumvent, iii) all claimed computations must actually be run.  
    - I switched to GPT as the main driver of the project, relegating Gemini agents to code audits, leading to consistent higher-quality and more scientifically rigourous outputs
    - Initial efforts of the GPT model towards an analytical proof (aided by computational checks) were not unpromising, but led to a rapid proliferation of code, dead-end directions, and increasingly baroque status reports. This was in part due to chasing transient patterns in low-level data. 
    - I realised that the model (and myself) was placing too much emphasis on the status of the protected baseline level 10 results, causing us to initially overlook and later fail to properly appreciate the fact that beyond level 5 ambiguities appear leading to alternative solutions.

Faced with the complexity of the proof effort here, and feeling that some key features of the initial computational evidence had been overlooked, I decided that it would be useful to start a competing investigation using an alternative model.
    
### Phase 3. A complementary restart 

- **Setup:** new project repository with protected files and agent instructions/skills as before.
    - I decided to work solely with Claude Opus 5.5 in Claude Code in order to test the capabilities of this line of models and gain an independent analysis
    - Initially I did not reveal any previous results to the model other than the existence of the baseline results up to level 10, which I emphasised may not be unique
- **Complementary progress:**
    - Earlier computational results were quickly reproduced and extended to higher levels: ultimately up to level 40
    - I checked that the model's initial results and proposed proof directions were in line with what had appeared in previous phases
    - I introduced a write-up (PREVIOUS_EFFORTS.md) with key comments about previously-tried methods, warning against pitfalls encountered, allowing the model to refine its approach
    - The model investigated some possible avenues for analytic progress, the efforts so far appear somewhat cleaner than the approach of the GPT model in phase 2, however I have not pushed this program very far nor verified the methods
    
## Usage level

These investigations have been carried out using mid-tier subscriptions (Google AI Pro, ChatGPT Plus, Claude Pro).
I have only been working intermittently on the problem (5-10 hours per week over 7-8 weeks) due to other committments.
Nonetheless I have frequently found myself hitting usage limits in the middle of reasoning efforts, leading to frustrating stop-start progress.

## Lessons for future research

1. Freeze a human-verified baseline early, and make it read-only to agents.
2. Write integrity rules down, and audit against them.
3. Restart cleanly with a written summary of lessons, rather than continuing a contaminated thread.
4. Use a second model for audit.
5. Watch for refinement spirals, and judge when to steer and when to prune.
6. Increase the subscription level and budget for additional compute.
