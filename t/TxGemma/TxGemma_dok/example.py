import os
import json
from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig, pipeline
from huggingface_hub import hf_hub_download, login, HfFolder

# Save the token under cache before login to Hugging Face
    try:
        token = HfFolder.get_token()
    except:
        print("Warning: Could not load cached token. Please run 'huggingface-cli login' or set HF_TOKEN environment variable.")

def run_example(MODEL_VARIANT, METHOD):
    print("This is an example of using the TxGemma model with Hugging Face Transformers.")

    # 1. Load model from Hugging Face Hub
    model_id = f"google/txgemma-{MODEL_VARIANT}"

    if MODEL_VARIANT == "2b-predict":
        additional_args = {}
    else:
        additional_args = {
            "quantization_config": BitsAndBytesConfig(load_in_8bit=True)
        }

    if METHOD == "direct":
        model = AutoModelForCausalLM.from_pretrained(
            model_id,
            device_map = "auto",
            **additional_args,
        )
        tokenizer = AutoTokenizer.from_pretrained(model_id)

    elif METHOD == "pipeline":
        pipe = pipeline(
        "text-generation",
        model = model,
        tokenizer = tokenizer,
        )

    # 2. Format prompts for therapeutic tasks
    # 2.1 Load prompt template
    tdc_prompts_filepath = hf_hub_download(
        repo_id=model_id,
        filename="tdc_prompts.json",
    )

    with open(tdc_prompts_filepath, "r") as f:
        tdc_prompts_json = json.load(f)

    # 2.2 Prepare sample prompt
    # Set example task and input
    task_name = "BBB_Martins"
    input_type = "{Drug SMILES}"
    drug_smiles = "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21"

    TDC_PROMPT = tdc_prompts_json[task_name].replace(input_type, drug_smiles)
    print("Formatted prompt:\n")
    print(TDC_PROMPT)

    # 3. Explore conversational capabilities with TxGemma-Chat
    # Gemini chat template: https://ai.google.dev/gemma/docs/core/prompt-structure?hl=zh-cn

    #3.1 Ask questions in a multi-turn conversation
    if METHOD == "direct":
        questions = [
            TDC_PROMPT,  # Initial question is a predictive task from TDC
            "Explain your reasoning based on the molecule structure."
        ]

        messages = []

        print("\n\n" + "="*50 + "\n")
        for question in questions:
            print(f"USER:\n{question}\n\n" + "-"*50 + "\n")
            messages.append(
                { "role": "user", "content": question },
            )
            # Apply the tokenizer's built-in chat template
            inputs = tokenizer.apply_chat_template(messages, tokenize=True, add_generation_prompt=True, return_tensors="pt")
            outputs = model.generate(input_ids=inputs.to("cuda"), max_new_tokens=512)
            response = tokenizer.decode(outputs[0, len(inputs[0]):], skip_special_tokens=True)
            print(f"TXGEMMA:\n{response}\n\n" + "="*50 + "\n")
            messages.append(
                { "role": "assistant", "content": response},
            )

    elif METHOD == "pipeline":
        questions = [
            TDC_PROMPT,  # Initial question is a predictive task from TDC
            "Explain your reasoning based on the molecule structure."
        ]

        messages = []

        print("\n\n" + "="*50 + "\n")
        for question in questions:
            print(f"USER:\n{question}\n\n" + "-"*50 + "\n")
            messages.append(
                { "role": "user", "content": question },
            )
            outputs = pipe(messages, max_new_tokens=512)
            messages = outputs[0]["generated_text"]
            response = messages[-1]["content"].strip()
            print(f"TXGEMMA:\n{response}\n\n" + "="*50 + "\n")



if __name__ == "__main__":
    MODEL_VARIANT = "9b-chat"  # @param ["2b-predict", "9b-chat", "9b-predict", "27b-chat", "27b-predict"]
    METHOD = "direct"  # @param ["direct", "pipeline"]
    run_example(MODEL_VARIANT, METHOD)
    print("Example completed successfully.")